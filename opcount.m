function varargout = opcount(action,varargin)

% -----------------------------------------------------
% -- OPCOUNT: MATLAB Arithmetic Operation Counting and Dataflow Analysis Tool
% -- Author: MD Salman Sakib
% -- Email: md_salman.sakib@okstate.edu

% -----------------------------------------------------
% -- Collaboration: Claude (AI assistant)
% -- Research paper: Energy Efficient Near-Threshold Voltage Massive
%                  MIMO Computing
% -----------------------------------------------------

% This tool is designed to analyze the computational complexity of MATLAB functions by dynamically tracing and counting arithmetic operations. 
% It provides detailed operation statistics and dataflow graph visualizations.

% -----------------------------------------------------
% The tool supports the following features:
% (1) Arithmetic operation counting: real and complex
%     multiplication, addition, division, and square root.
% (2) Dynamic function tracing: execution of numerical
%     algorithms, loops, branches, and nested function calls.
% (3) Dataflow analysis: visualization of arithmetic
%     dependencies, operation counts, and execution paths.
% (4) Computational complexity analysis: evaluation of
%     numerical and massive MIMO detection algorithms.


% -----------------------------------------------------
% Please cite the following paper if you use this tool:
% MD Salman Sakib, Mehdi Safarpour, Pejman Ghasemzadeh, Shahriar Shahabddin, "Energy Efficient Near-Threshold Voltage Massive MIMO Computing"
%
% -----------------------------------------------------

% Purpose: Main user-facing controller that handles tic, toc, graph, and reset,
% replays the enclosed function, counts arithmetic, and produces statistics.
% Inputs: action selects tic, toc, graph, or reset; varargin provides optional
% graph arguments.
% Output(s): varargout contains statistics, graph data, or an empty result.
%
% Maintainer note: Explanatory comments are included throughout the code.
% Executable statements are preserved from the original implementation.
%OPCOUNT General dynamic arithmetic counter for ordinary MATLAB functions.
%
% Usage:
%   opcount tic
%   [out1,out2] = my_function(in1,in2);
%   stats = opcount('toc');
%   opcount graph
%
% The enclosed user function is replayed automatically with traced versions
% of its actual inputs. Arithmetic executed in the selected function and in
% user-defined functions called from it is counted without modifying the
% user's expressions. The replay follows the actual loop counts and branch
% decisions produced by those inputs.
%
% Counted operation classes:
%   CMUL, CADD, RMUL, RADD, DIV, SQRT.
%
% The graph is intentionally coarse and algorithm-independent. Runtime
% operations are grouped by source statement, operator position, call path,
% and loop context. Repeated loop executions collapse into one logical loop
% body. Elementwise chains on one source statement may be merged into a
% compound node. Each node shows source code, line number, operation class,
% runtime dimensions, total operation count, and execution count.
%
% The region must contain exactly one top-level function call. Nested and
% local user-defined calls made by that function are traced recursively.
%
% Important: opcount('toc') replays the selected numerical function once.
% RNG state is restored, but other external side effects should not be put
% inside the selected function.
%
% Runtime safety: the runtime no longer places wrappers named zeros,
% ones, false, randn, etc. on the MATLAB path.
%
% Local-function replay: generated factory files are compact. Source
% locations are translated through an explicit generated-line to original-line
% map, so a local function may begin anywhere in the user's file.
%
% V28 preserves traced indexed assignments into zero-initialized arrays:
%   - bare zeros(...) constructors are promoted only in the temporary replay copy;
%   - the user's original source and original execution are unchanged;
%   - indexed updates retain dependency/support metadata across loop iterations;
%   - no arithmetic cost model from V27 is changed.
%
% Decomposition-aware structural tracing remains enabled:
%   - Cholesky, classical Gram-Schmidt QR, LU, and LDL operation models;
%   - output-specific lower/upper/diagonal/permutation support metadata;
%   - transpose propagates lower <-> upper structure into triangular solves;
%   - helper operations inherit the caller loop without creating duplicate regions;
%   - decomposition masks are shape-normalized and CHOL uses guaranteed triangular support.
% exact zero structure is inferred for numeric inputs and then propagated through
% the trace, so an upper/lower-triangular argument remains triangular even when
% it was formed outside the counted function. Triangular inversion follows the
% supplied recurrence directly: diagonal reciprocals are computed once, each
% off-diagonal inner product skips guaranteed zeros, and the result is scaled by
% the stored reciprocal diagonal. Counts remain literal to the source; repeated
% loop-invariant work is counted rather than silently hoisted.
%
% MATLAB syntax note: use stats = opcount('toc'). Command syntax cannot
% assign the output of "opcount toc".

% Remove any path-shadowing constructor wrappers left by the failed v9
% prototype before this function uses logical constants or array helpers.
purge_legacy_runtime();

persistent state last_flow
if nargout > 1
    error('opcount:TooManyOutputs','opcount returns at most one output.');
end
if isempty(state), state = empty_state(); end
if isempty(last_flow), last_flow = []; end
if nargin == 0, action = 'toc'; end
if isstring(action)
    if ~isscalar(action), error('opcount:InvalidAction','Action must be scalar.'); end
    action = char(action);
end
if ~ischar(action)
    error('opcount:InvalidAction','Use opcount tic, opcount(''toc''), opcount graph, or opcount reset.');
end
action = lower(strtrim(action));

switch action
    case 'tic'
        if state.active
            error('opcount:AlreadyActive','An opcount region is already active.');
        end
        caller = get_caller();
        state = empty_state();
        state.active = true;
        state.caller_file = caller.file;
        state.caller_name = caller.name;
        state.start_line = caller.line;
        state.timer = tic;
        try
            state.rng_state = rng;
        catch
            state.rng_state = [];
        end
        result = [];

    case 'toc'
        if ~state.active
            error('opcount:NotActive','No active region exists. Start with "opcount tic".');
        end
        caller = get_caller();
        snapshot = state;
        state = empty_state();
        if ~same_file(caller.file,snapshot.caller_file) || ~strcmp(caller.name,snapshot.caller_name)
            error('opcount:CallerChanged','tic and toc must execute in the same MATLAB script or function.');
        end
        if caller.line <= snapshot.start_line
            error('opcount:InvalidRegion','The toc marker must occur after the tic marker.');
        end
        elapsed_seconds = toc(snapshot.timer);
        call = find_single_enclosed_call(snapshot.caller_file,snapshot.start_line,caller.line);
        inputs = cell(1,numel(call.input_expressions));
        for kk=1:numel(call.input_expressions)
            try
                inputs{kk}=evalin('caller',call.input_expressions{kk});
            catch ME
                error('opcount:InputEvaluationFailed', ...
                    'Could not evaluate input expression "%s": %s', ...
                    call.input_expressions{kk},ME.message);
            end
        end
        outputs = cell(1,numel(call.output_expressions));
        for kk=1:numel(call.output_expressions)
            expr=strtrim(call.output_expressions{kk});
            if isempty(expr)||strcmp(expr,'~')
                outputs{kk}=[];
            else
                try
                    outputs{kk}=evalin('caller',expr);
                catch ME
                    error('opcount:OutputEvaluationFailed', ...
                        'Could not evaluate output expression "%s": %s',expr,ME.message);
                end
            end
        end
        % Resolve a callable handle for the enclosed function. MATLAB's
        % evalin cannot construct a handle to a function defined in the caller.
        % For a local function, OPCOUNT therefore creates a temporary factory
        % file containing the caller's local-function section and obtains a true
        % scoped handle through localfunctions(). Standalone/package functions
        % are resolved normally with str2func.
        runtime_root = ensure_trace_class();
        try
            [function_handle,source_maps] = resolve_replay_handle(call,snapshot.caller_file,runtime_root);
        catch ME
            error('opcount:FunctionHandleCreationFailed', ...
                'Could not resolve a replay handle for "%s": %s', ...
                call.function_name,ME.message);
        end
        metadata = call;
        metadata.caller_file = snapshot.caller_file;
        metadata.caller_name = snapshot.caller_name;
        metadata.elapsed_seconds = elapsed_seconds;
        metadata.source_region = [snapshot.start_line+1 caller.line-1];

        [raw_nodes,replay_outputs] = replay_with_tracing(call,function_handle,inputs,snapshot.rng_state,runtime_root,source_maps); %#ok<ASGLU>
        [stats,last_flow] = build_results(raw_nodes,call,inputs,outputs,metadata);
        print_report(stats);
        result = stats;

    case 'graph'
        if isempty(last_flow)
            error('opcount:NoGraph','No graph is available. Complete a successful tic/toc region first.');
        end
        plot_dataflow(last_flow,varargin{:});
        result = last_flow;

    case 'reset'
        state = empty_state();
        last_flow = [];
        clear_runtime();
        result = [];

    otherwise
        error('opcount:UnknownAction','Unknown action "%s".',action);
end

if nargout > 0, varargout{1} = result; end
end

function state = empty_state()
state = struct('active',false,'caller_file','','caller_name','', ...
    'start_line',[],'timer',[],'rng_state',[]);
end

function caller = get_caller()
stack = dbstack('-completenames');
if isempty(stack), error('opcount:NoCaller','Call opcount from a saved MATLAB file.'); end
internal = stack(1).file;
caller = [];
for k=2:numel(stack)
    if ~same_file(stack(k).file,internal)
        caller = stack(k); break
    end
end
if isempty(caller) || isempty(caller.file)
    error('opcount:UnsavedCaller','The calling script or function must be saved to a file.');
end
end

function purge_legacy_runtime()
% v9 added general constructor names (false, zeros, randn, ...) to the
% beginning of the MATLAB path. An out-of-memory interruption could leave
% that directory active and poison later MATLAB calls. Remove it before any
% OPCOUNT state is initialized. The current runtime contains only an @class
% directory and never shadows MATLAB constructors.
legacy_roots = {fullfile(tempdir,'opcount_runtime_20260804a')};
parts = strsplit(path,pathsep);
for k = 1:numel(legacy_roots)
    if any(strcmp(parts,legacy_roots{k}))
        try, rmpath(legacy_roots{k}); catch, end
    end
end
global OPCOUNT_RUNTIME_V9
OPCOUNT_RUNTIME_V9 = [];
end

function tf = same_file(a,b)
tf = ~isempty(a) && ~isempty(b) && strcmpi(char(a),char(b));
end

function values = evaluate_expressions(expressions,workspace,kind)
values = cell(1,numel(expressions));
for k=1:numel(expressions)
    try
        values{k} = evalin(workspace,expressions{k});
    catch ME
        error('opcount:ExpressionEvaluationFailed','Could not evaluate %s expression "%s": %s',kind,expressions{k},ME.message);
    end
end
end

function outputs = evaluate_outputs(expressions)
outputs = cell(1,numel(expressions));
for k=1:numel(expressions)
    expr = strtrim(expressions{k});
    if isempty(expr) || strcmp(expr,'~')
        outputs{k} = [];
    else
        try
            outputs{k} = evalin('caller',expr);
        catch ME
            error('opcount:OutputEvaluationFailed','Could not evaluate output "%s": %s',expr,ME.message);
        end
    end
end
end

function [fh,source_maps] = resolve_replay_handle(call,caller_file,runtime_root)
name = char(call.function_name);
source_maps = empty_source_maps();
info = inspect_caller_functions(caller_file,name);

if strcmp(info.kind,'local')
    [factory_name,source_map] = build_local_function_factory(caller_file,runtime_root,info.local_section_start);
    source_maps = source_map;
    factory = str2func(factory_name);
    fh = factory(name);
    return
elseif strcmp(info.kind,'nested')
    error('opcount:NestedTargetUnsupported', ...
        ['The enclosed target "%s" is a nested function. A nested target can ' ...
         'capture variables from its parent workspace and cannot be reconstructed ' ...
         'safely from the source file. Enclose its parent function or move the ' ...
         'nested target to a local/standalone function.'],name);
end

% Main, standalone, package, and class-qualified functions.
fh = str2func(name);
finfo = functions(fh);
if isempty(which(name)) && (~isfield(finfo,'file') || isempty(finfo.file))
    error('opcount:FunctionNotFound','Function "%s" is not visible on the MATLAB path.',name);
end
end

function info = inspect_caller_functions(file_name,target_name)
info = struct('kind','external','local_section_start',0,'target_line',0,'main_end',0);
if isempty(file_name) || exist(file_name,'file') ~= 2, return, end
lines = regexp(fileread(file_name),'\r\n|\n|\r','split');
if isempty(lines), return, end

[first_code,first_word] = first_code_token(lines);
is_function_file = first_code > 0 && strcmp(first_word,'function');
main_end = 0;
if is_function_file
    main_end = matching_main_function_end(lines,first_code);
    if main_end == 0
        error('opcount:SourceParseFailed','Could not locate the end of the main function in %s.',file_name);
    end
end

[target_line,all_function_lines] = find_function_declaration(lines,target_name);
if target_line == 0, return, end
info.target_line = target_line;
info.main_end = main_end;

if is_function_file
    if target_line == first_code
        info.kind = 'main';
    elseif target_line > main_end
        info.kind = 'local';
        info.local_section_start = main_end + 1;
    else
        info.kind = 'nested';
    end
else
    % In a script, all function definitions are local to the script. Preserve
    % the complete local-function section beginning at the first declaration.
    info.kind = 'local';
    info.local_section_start = min(all_function_lines);
end
end

function [line_number,word] = first_code_token(lines)
line_number = 0; word = '';
in_block_comment = false;
for k = 1:numel(lines)
    raw = strtrim(lines{k});
    if in_block_comment
        if startsWith(raw,'%}'), in_block_comment = false; end
        continue
    end
    if startsWith(raw,'%{'), in_block_comment = true; continue, end
    code = strtrim(strip_comment(lines{k}));
    if isempty(code), continue, end
    token = regexp(code,'^([A-Za-z]\w*)','tokens','once');
    if isempty(token), continue, end
    line_number = k; word = lower(token{1}); return
end
end

function main_end = matching_main_function_end(lines,start_line)
% Locate the explicit end of the primary function. MATLAB permits compact
% control blocks on one physical line, for example:
%   if condition, value = 1; end
% Counting only the first token on that line leaves the if unmatched and can
% make every later local function appear nested. Parse every top-level
% comma/semicolon-delimited control statement instead.
main_end = 0; depth = 0; in_block_comment = false;
openers = {'function','if','for','parfor','while','switch','try','spmd', ...
    'classdef','properties','methods','events','enumeration','arguments'};
for k = start_line:numel(lines)
    raw = strtrim(lines{k});
    if in_block_comment
        if startsWith(raw,'%}'), in_block_comment = false; end
        continue
    end
    if startsWith(raw,'%{'), in_block_comment = true; continue, end
    code = strip_comment(lines{k});
    if isempty(strtrim(code)), continue, end
    statements = split_top_level_control_statements(code);
    for q = 1:numel(statements)
        statement = strtrim(statements{q});
        token = regexp(statement,'^([A-Za-z]\w*)','tokens','once');
        if isempty(token), continue, end
        word = lower(token{1});
        if any(strcmp(word,openers))
            depth = depth + 1;
        elseif strcmp(word,'end')
            depth = depth - 1;
            if depth == 0
                main_end = k; return
            elseif depth < 0
                return
            end
        end
    end
end
end

function [target_line,function_lines] = find_function_declaration(lines,target_name)
target_line = 0; function_lines = [];
in_block_comment = false;
for k = 1:numel(lines)
    raw = strtrim(lines{k});
    if in_block_comment
        if startsWith(raw,'%}'), in_block_comment = false; end
        continue
    end
    if startsWith(raw,'%{'), in_block_comment = true; continue, end
    code = strtrim(strip_comment(lines{k}));
    if isempty(regexp(code,'^function(?:\s|$)','once')), continue, end
    function_lines(end+1) = k; %#ok<AGROW>
    declared = function_name_from_header(code);
    if strcmp(declared,target_name), target_line = k; end
end
end

function name = function_name_from_header(code)
name = '';
rest = strtrim(regexprep(code,'^function(?:\s|$)',''));
[eq,~] = top_level_assignment(rest);
if eq > 0, rest = strtrim(rest(eq+1:end)); end
match = regexp(rest,'^([A-Za-z]\w*)','tokens','once');
if ~isempty(match), name = match{1}; end
end

function maps = empty_source_maps()
maps = struct('generated_file',{},'original_file',{}, ...
    'generated_start',{},'original_start',{});
end

function name = local_factory_name(source_file)
[~,base] = fileparts(source_file);
base = regexprep(base,'[^A-Za-z0-9_]','_');
if isempty(base) || ~isletter(base(1)), base = ['f_' base]; end
if numel(base) > 20, base = base(1:20); end
bytes = uint8(lower(char(source_file)));
h = uint32(2166136261);
for k = 1:numel(bytes)
    h = bitxor(h,uint32(bytes(k)));
    h = uint32(mod(uint64(h)*uint64(16777619),uint64(4294967296)));
end
name = sprintf('opcount_factory_%s_%08x',base,h);
end

function [factory_name,source_map] = build_local_function_factory(source_file,runtime_root,section_start)
factory_name = local_factory_name(source_file);
factory_file = fullfile(runtime_root,[factory_name '.m']);
lines = regexp(fileread(source_file),'\r\n|\n|\r','split');
if section_start < 1 || section_start > numel(lines)
    error('opcount:LocalSectionNotFound','Could not locate the local-function section in %s.',source_file);
end

% Keep the generated factory compact. Absolute source-line preservation is
% deliberately not attempted; the tracing runtime translates generated line
% numbers back to the original file through source_map.
preamble = {
    sprintf('function fh = %s(targetName)',factory_name)
    'fcns = localfunctions;'
    'fh = [];'
    'for opcountFactoryIndex = 1:numel(fcns)'
    '    opcountFactoryName = func2str(fcns{opcountFactoryIndex});'
    '    opcountFactoryParts = regexp(opcountFactoryName,''[>/]'',''split'');'
    '    if strcmp(opcountFactoryParts{end},targetName)'
    '        fh = fcns{opcountFactoryIndex};'
    '        return'
    '    end'
    'end'
    'error(''opcount:LocalFunctionNotFound'',''Local function "%s" was not found in the generated factory.'',targetName);'
    'end'
    };

generated_start = numel(preamble) + 1;
replay_lines = lines(section_start:end).';
% Numeric zero arrays created with ordinary numeric dimensions would normally
% remain MATLAB doubles during replay.  If traced values are later assigned
% into such arrays, MATLAB converts the values back to doubles and loses the
% arithmetic/dependency metadata.  Rewrite only bare ZEROS(...) constructors
% in the temporary replay copy so indexed assignments remain trace-aware.
% The user's original source is never modified, and line numbers are preserved.
replay_lines = rewrite_replay_zero_initializers(replay_lines,trace_class_name());
out = [preamble; replay_lines];
source = strjoin(out,newline);
source_map = struct('generated_file',factory_file,'original_file',source_file, ...
    'generated_start',generated_start,'original_start',section_start);

write_file = true;
if exist(factory_file,'file') == 2
    try, write_file = ~strcmp(fileread(factory_file),source); catch, end
end
if write_file
    fid = fopen(factory_file,'w');
    if fid < 0, error('opcount:FactoryWriteFailed','Could not create %s.',factory_file); end
    cleaner = onCleanup(@()fclose(fid)); %#ok<NASGU>
    fwrite(fid,source,'char');
    clear cleaner
    rehash
end
end

function lines = rewrite_replay_zero_initializers(lines,class_name)
% Preserve trace metadata for arrays used as indexed assignment containers.
replacement = [class_name '.trace_zeros'];
for k = 1:numel(lines)
    code = strtrim(strip_comment(lines{k}));
    % Do not rewrite a user function declaration named ZEROS.
    if ~isempty(regexp(code,'^function(?:\s|$).*\bzeros\s*\(','once'))
        continue
    end
    lines{k} = replace_bare_function_call(lines{k},'zeros',replacement);
end
end

function out = replace_bare_function_call(line,name,replacement)
% Replace NAME(...) only when NAME is a bare function call outside strings
% and comments.  Object methods such as obj.zeros(...) are left untouched.
out = '';
in_string = false;
i = 1;
name_len = numel(name);
while i <= numel(line)
    ch = line(i);
    if ch == ''''
        if in_string && i < numel(line) && line(i+1) == ''''
            out = [out ch line(i+1)]; %#ok<AGROW>
            i = i + 2;
            continue
        end
        prev = previous_nonspace(out);
        if ~in_string && ~isempty(prev) && ...
                (isstrprop(prev,'alphanum') || any(prev == '_])}'))
            % Conjugate-transpose operator, not the start of a string.
            out = [out ch]; %#ok<AGROW>
            i = i + 1;
            continue
        end
        in_string = ~in_string;
        out = [out ch]; %#ok<AGROW>
        i = i + 1;
        continue
    end
    if ch == '%' && ~in_string
        out = [out line(i:end)]; %#ok<AGROW>
        break
    end
    if ~in_string && i + name_len - 1 <= numel(line) && ...
            strcmp(line(i:i+name_len-1),name)
        before_ok = (i == 1) || ...
            ~(isstrprop(line(i-1),'alphanum') || line(i-1) == '_' || line(i-1) == '.');
        after_index = i + name_len;
        after_ok = after_index > numel(line) || ...
            ~(isstrprop(line(after_index),'alphanum') || line(after_index) == '_');
        probe = after_index;
        while probe <= numel(line) && isspace(line(probe)), probe = probe + 1; end
        call_ok = probe <= numel(line) && line(probe) == '(';
        if before_ok && after_ok && call_ok
            out = [out replacement]; %#ok<AGROW>
            i = after_index;
            continue
        end
    end
    out = [out ch]; %#ok<AGROW>
    i = i + 1;
end
end

function [raw_nodes,replay_outputs] = replay_with_tracing(call,fh,inputs,rng_state,runtime_root,source_maps)
global OPCOUNT_RUNTIME_V28
OPCOUNT_RUNTIME_V28 = struct();
OPCOUNT_RUNTIME_V28.active = true;
OPCOUNT_RUNTIME_V28.next_id = 1;
OPCOUNT_RUNTIME_V28.nodes = empty_raw_nodes();
OPCOUNT_RUNTIME_V28.source_cache = containers.Map('KeyType','char','ValueType','any');
OPCOUNT_RUNTIME_V28.source_maps = source_maps;
cleanup_path = onCleanup(@()remove_runtime_path(runtime_root)); %#ok<NASGU>

wrapped = cell(size(inputs));
for k=1:numel(inputs)
    wrapped{k} = wrap_value(inputs{k});
end

current_rng = [];
try, current_rng = rng; catch, end
cleanup_rng = onCleanup(@()restore_rng(current_rng)); %#ok<NASGU>
if ~isempty(rng_state)
    try, rng(rng_state); catch, end
end

nout = numel(call.output_expressions);
replay_outputs = cell(1,nout);
try
    if nout == 0
        fh(wrapped{:});
    else
        [replay_outputs{:}] = fh(wrapped{:});
    end
catch ME
    raw_nodes = OPCOUNT_RUNTIME_V28.nodes;
    detail = ME.message;
    if ~isempty(ME.stack)
        frame = ME.stack(1);
        detail = sprintf('%s [at %s, line %d]',detail,frame.name,frame.line);
    end
    clear_runtime();
    error('opcount:ReplayFailed', ...
        ['Automatic tracing failed while replaying "%s". The current general analyzer supports ordinary numerical MATLAB code, nested user functions, loops, branches, indexing, common matrix operations, and common linear-algebra built-ins. Original error: %s'], ...
        call.function_name,detail);
end
raw_nodes = OPCOUNT_RUNTIME_V28.nodes;
OPCOUNT_RUNTIME_V28.active = false;
end

function restore_rng(s)
if ~isempty(s), try, rng(s); catch, end, end
end

function clear_runtime()
global OPCOUNT_RUNTIME_V28
OPCOUNT_RUNTIME_V28 = [];
end

function wrapped = wrap_value(value)
class_name = trace_class_name();
if isnumeric(value) || islogical(value)
    % Infer exact support for top-level numeric inputs. This lets matrices
    % constructed outside the counted function (for example with TRIU/TRIL)
    % retain their diagonal/triangular structure during replay.
    support = (value ~= 0);
    unit_mask = isreal(value) & ((value == 1) | (value == -1));
    wrapped = feval(class_name,value,builtin('zeros',1,0),false,support,unit_mask);
elseif isstruct(value)
    wrapped = value;
    fields = fieldnames(value);
    for idx=1:numel(value)
        for k=1:numel(fields)
            wrapped(idx).(fields{k}) = wrap_value(value(idx).(fields{k}));
        end
    end
elseif iscell(value)
    wrapped = value;
    for k=1:numel(value), wrapped{k}=wrap_value(value{k}); end
else
    wrapped = value;
end
end

function name = trace_class_name()
name = 'opcount_trace_20260805j';
end

function root = ensure_trace_class()
name = trace_class_name();
root = fullfile(tempdir,'opcount_runtime_20260805j');
class_dir = fullfile(root,['@' name]);
class_file = fullfile(class_dir,[name '.m']);
if exist(class_dir,'dir') ~= 7, mkdir(class_dir); end
source = embedded_trace_source();
write_file = true;
if exist(class_file,'file') == 2
    try, write_file = ~strcmp(fileread(class_file),source); catch, end
end
if write_file
    fid=fopen(class_file,'w');
    if fid<0, error('opcount:TraceClassWriteFailed','Could not create the temporary tracing class.'); end
    cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fwrite(fid,source,'char');
    clear cleaner
end
if isempty(strfind([path pathsep],[root pathsep])) %#ok<STREMP>
    old_warning=warning('off','MATLAB:dispatcher:nameConflict');
    cleanup_warning=onCleanup(@()warning(old_warning)); %#ok<NASGU>
    addpath(root,'-begin');
end
if exist(name,'class') ~= 8
    rehash;
end
end

function remove_runtime_path(root)
if isempty(root), return, end
parts=strsplit(path,pathsep);
if any(strcmp(parts,root))
    try, rmpath(root); catch, end
end
end

function source = embedded_trace_source()
stack = dbstack('-completenames');
this_file = stack(1).file;
if isempty(this_file), error('opcount:SelfLocationFailed','Could not locate opcount.m.'); end
text = fileread(this_file);
begin_marker = '%<OPCOUNT_TRACE_CLASS_BEGIN>';
end_marker = '%<OPCOUNT_TRACE_CLASS_END>';
a = strfind(text,begin_marker); b = strfind(text,end_marker);
if isempty(a) || isempty(b) || b(end)<=a(end)
    error('opcount:EmbeddedClassMissing','The embedded tracing class is missing from opcount.m.');
end
block = text(a(end)+numel(begin_marker):b(end)-1);
lines = regexp(block,'\r\n|\n|\r','split');
for k=1:numel(lines)
    if startsWith(lines{k},'%|'), lines{k}=lines{k}(3:end); end
end
source = strjoin(lines,newline);
source = regexprep(source,'^\s*\n','');
end

function nodes = empty_raw_nodes()
template = struct('id',0,'parents',builtin('zeros',1,0),'counts',zero_counts(), ...
    'file','','function_name','','line',0,'expression','','operator','', ...
    'stack_frames',{{}},'input_dimensions',{{}},'output_dimensions','');
nodes = repmat(template,1,0);
end

function c = zero_counts()
c = struct('complex_mult',0,'complex_add',0,'real_mult',0,'real_add',0,'divisions',0,'square_roots',0, 'sort_comparisons',0);
end

function [stats,flow] = build_results(raw_nodes,call,inputs,outputs,metadata)
[groups,edges,feedback,loop_regions] = coarsen_nodes(raw_nodes,call.function_name);
stats = zero_counts();
for k=1:numel(raw_nodes), stats = add_counts(stats,raw_nodes(k).counts); end
stats.real_equivalent_mult = stats.real_mult + 4*stats.complex_mult;
stats.real_equivalent_add = stats.real_add + 2*stats.complex_add + 2*stats.complex_mult;
stats.function_name = call.function_name;
stats.source_region = metadata.source_region;
stats.elapsed_seconds = metadata.elapsed_seconds;
stats.input_sizes = cellfun(@size,inputs,'UniformOutput',false);
stats.output_sizes = cellfun(@size,outputs,'UniformOutput',false);
stats.call_statement = call.statement;

[depth,max_width,critical] = graph_metrics(groups,edges);
flow = struct();
flow.function_name = call.function_name;
flow.call_statement = call.statement;
flow.nodes = groups;
flow.edges = edges;
flow.feedback_edges = feedback;
flow.loop_regions = loop_regions;
flow.source_lines = build_source_line_summary(groups);
flow.metrics = struct('operation_nodes',numel(groups), ...
    'max_parallel_width',max_width,'critical_path_stages',critical);
flow.total_counts = stats;
flow.source_region = metadata.source_region;
stats.dataflow = flow.metrics;
end

function total = add_counts(total,c)
fields=fieldnames(total);
for k=1:numel(fields), total.(fields{k})=total.(fields{k})+c.(fields{k}); end
end

function [groups,edges,feedback,loop_regions] = coarsen_nodes(raw,root_function)
% Algorithm-independent coarsening:
%   1) identify a structural operator slot within each source statement;
%   2) merge the same slot across repeated loop executions, even when the
%      first pass changes real/complex type;
%   3) merge safe elementwise chains on the same source statement;
%   4) retain dependencies and explicit feedback edges.
if isempty(raw)
    groups=empty_group_nodes(); edges=builtin('zeros',0,2);
    feedback=builtin('zeros',0,2); loop_regions=empty_loop_regions(); return
end
[ordinals,slot_counts,loop_info,execution_weights] = operation_ordinals(raw,root_function);
keys=cell(1,numel(raw)); map=builtin('zeros',1,numel(raw));
groups=empty_group_nodes(); unique_keys={};
for k=1:numel(raw)
    statement_key=statement_identity(raw(k));
    call_key=call_path_key(raw(k).stack_frames,root_function);
    keys{k}=sprintf('%s|CALL:%s|SLOT:%d',statement_key,call_key,ordinals(k));
    idx=find(strcmp(unique_keys,keys{k}),1);
    if isempty(idx)
        unique_keys{end+1}=keys{k}; %#ok<AGROW>
        idx=numel(unique_keys);
        g=empty_group_node();
        g.id=idx; g.first_raw_id=raw(k).id; g.last_raw_id=raw(k).id;
        g.file=raw(k).file;
        g.function_name=raw(k).function_name; g.line=raw(k).line;
        g.expression=raw(k).expression; g.statement_key=statement_key;
        g.call_path=call_key; g.counts=zero_counts();
        g.slot_start=ordinals(k); g.slot_end=ordinals(k);
        g.statement_slots=slot_counts(k);
        if isempty(loop_info(k).id)
            g.loop_id='';
        else
            % The active source loop is the region identity. A helper called
            % from that loop inherits the same region; CALL path still remains
            % in the node key, so helper and caller operations stay distinct.
            g.loop_id=loop_info(k).id;
        end
        g.loop_label=loop_info(k).label;
        g.loop_line=loop_info(k).line;
        g.loop_end_line=loop_info(k).end_line;
        groups(end+1)=g; %#ok<AGROW>
    end
    map(k)=idx;
    groups(idx).executions=groups(idx).executions+execution_weights(k);
    groups(idx).last_raw_id=max(groups(idx).last_raw_id,raw(k).id);
    groups(idx).counts=add_counts(groups(idx).counts,raw(k).counts);
    groups(idx).operators=append_unique(groups(idx).operators,normalized_operator(raw(k)));
    input_text=strjoin(raw(k).input_dimensions,' x ');
    groups(idx).input_variants=append_unique(groups(idx).input_variants,input_text);
    groups(idx).output_variants=append_unique(groups(idx).output_variants,raw(k).output_dimensions);
end
for k=1:numel(groups), groups(k)=finalize_group(groups(k)); end

all_edges=builtin('zeros',0,2);
for k=1:numel(raw)
    child=map(k);
    for p=raw(k).parents
        rp=find([raw.id]==p,1);
        if isempty(rp), continue, end
        parent=map(rp);
        all_edges(end+1,:)=[parent child]; %#ok<AGROW>
    end
end
if isempty(all_edges)
    edges=builtin('zeros',0,2); feedback=builtin('zeros',0,2);
else
    all_edges=unique(all_edges,'rows','stable');
    edges=builtin('zeros',0,2); feedback=builtin('zeros',0,2);
    for k=1:size(all_edges,1)
        a=all_edges(k,1); b=all_edges(k,2);
        if a==b || groups(a).first_raw_id>=groups(b).first_raw_id
            feedback(end+1,:)=[a b]; %#ok<AGROW>
        else
            edges(end+1,:)=[a b]; %#ok<AGROW>
        end
    end
end
[groups,edges,feedback]=merge_elementwise_chains(groups,edges,feedback);
loop_regions=collect_loop_regions(groups);
edges=add_loop_exit_barriers(groups,edges,loop_regions);
end

function [ordinals,slot_counts,loop_info,execution_weights] = operation_ordinals(raw,root_function)
% Find the repeating operator sequence for each source statement. Operator
% identity, rather than real/complex counts, defines the period. This lets
% the first loop pass differ in numeric type without creating extra nodes.
N=numel(raw); ordinals=ones(1,N); slot_counts=ones(1,N); execution_weights=ones(1,N);
loop_info=repmat(empty_loop_info(),1,N);
base_keys=cell(1,N);
for k=1:N
    base_keys{k}=sprintf('%s|CALL:%s',statement_identity(raw(k)),call_path_key(raw(k).stack_frames,root_function));
    loop_info(k)=nearest_loop_context(raw(k).stack_frames,root_function);
end
unique_lines=unique(base_keys,'stable');
for u=1:numel(unique_lines)
    idx=find(strcmp(base_keys,unique_lines{u}));
    L=numel(idx); ops=cell(1,L);
    for q=1:L, ops{q}=normalized_operator(raw(idx(q))); end
    p=L; found_period=false;
    in_loop=~isempty(loop_info(idx(1)).id);
    if in_loop
        for candidate=1:L
            if mod(L,candidate)~=0, continue, end
            repeats=L/candidate;
            if repeats<2, continue, end
            matches=true;
            for q=1:L
                if ~strcmp(ops{q},ops{mod(q-1,candidate)+1}), matches=false; break, end
            end
            if matches, p=candidate; found_period=true; break, end
        end
    end
    if in_loop && ~found_period && L>1
        % The operator sequence can vary by iteration (for example A^i).
        % Collapse the complete source statement into one coarse loop node.
        % Estimate the number of loop passes from runtime operator frequency
        % divided by the operator's syntactic occurrences in the statement.
        % This remains algorithm-independent and is robust when early passes
        % contain zero-cost identity/zero operations that are not graph nodes.
        p=1;
        pass_estimate=1;
        unique_ops=unique(ops,'stable');
        expr=raw(idx(1)).expression;
        for z=1:numel(unique_ops)
            runtime_count=sum(strcmp(ops,unique_ops{z}));
            static_count=static_operator_occurrences(expr,unique_ops{z});
            if static_count>0 && mod(runtime_count,static_count)==0
                pass_estimate=max(pass_estimate,runtime_count/static_count);
            end
        end
        execution_weights(idx)=0;
        execution_weights(idx(end))=pass_estimate;
    end
    for q=1:L
        ordinals(idx(q))=mod(q-1,p)+1;
        slot_counts(idx(q))=p;
    end
end
end


function n=static_operator_occurrences(expression,op)
% Conservative syntactic occurrence count used only to estimate loop passes
% when one statement executes a nonuniform operator sequence across passes.
text=strip_comment(char(expression));
op=lower(char(op));
switch op
    case {'inv','pinv','chol','qr','lu','norm','abs','sqrt','sum','mean','prod','dot'}
        n=numel(regexp(lower(text),['\<' op '\s*\('],'match'));
    case 'mpower'
        n=numel(regexp(text,'(?<!\.)\^(?!\^)','match'));
    case 'power'
        n=numel(regexp(text,'\.\^','match'));
    case 'mtimes'
        n=numel(regexp(text,'(?<!\.)\*(?!\*)','match'));
    case 'times'
        n=numel(regexp(text,'\.\*','match'));
    case {'rdivide','mrdivide'}
        n=numel(regexp(text,'(?<!\.)/(?!/)','match'));
    case 'plus'
        n=numel(regexp(text,'(?<![eE])\+','match'));
    case 'minus'
        n=numel(regexp(text,'(?<![eE])-','match'));
    otherwise
        n=0;
end
end

function key=statement_identity(n)
expr=regexprep(strtrim(n.expression),'\s+',' ');
key=sprintf('%s|%s|%d|%s',lower(n.file),lower(n.function_name),n.line,expr);
end

function key=call_path_key(frames,root_function)
% Exclude the innermost frame because file/function/line are already in the
% statement identity. Outer frames distinguish helper functions called from
% different source sites and preserve loop context.
key='';
frames=frames_within_trace_root(frames,root_function);
if isempty(frames)||~isstruct(frames)||numel(frames)<2, return, end
parts=cell(1,numel(frames)-1);
for k=2:numel(frames)
    parts{k-1}=sprintf('%s|%s|%d',lower(frames(k).file),lower(frames(k).function_name),frames(k).line);
end
key=strjoin(parts,'>');
end

function op=normalized_operator(n)
op='';
if isfield(n,'operator'), op=lower(strtrim(n.operator)); end
if isempty(op)
    c=n.counts;
    if c.divisions>0 && (c.complex_mult+c.real_mult)>0, op='solve';
    elseif c.square_roots>0 && (c.real_mult+c.real_add)>0, op='magnitude';
    elseif (c.complex_mult+c.real_mult)>0 && (c.complex_add+c.real_add)>0, op='mtimes';
    elseif (c.complex_mult+c.real_mult)>0, op='times';
    elseif (c.complex_add+c.real_add)>0, op='plus';
    elseif c.divisions>0, op='rdivide';
    elseif c.square_roots>0, op='sqrt';
    else, op='arithmetic'; end
end
end

function values=append_unique(values,value)
if isempty(value), return, end
if isempty(values)||~any(strcmp(values,value)), values{end+1}=value; end %#ok<AGROW>
end

function g=finalize_group(g)
g.operation_text=operation_text(g.counts);
g.count_text=count_text(g.counts);
g.operator_text=operator_sequence_text(g.operators);
g.dimension_text=group_dimension_text(g);
end

function [groups,edges,feedback]=merge_elementwise_chains(groups,edges,feedback)
% Merge only unambiguous, same-statement elementwise chains. Matrix
% products, inversions, and solves remain separate so parallel branches and
% critical-path structure stay visible. Example: minus -> abs -> power can
% become one compound distance node.
changed=true;
while changed
    changed=false;
    if isempty(groups)||isempty(edges), break, end
    for a=1:numel(groups)
        children=edges(edges(:,1)==a,2);
        if numel(children)~=1, continue, end
        b=children(1);
        parents_b=edges(edges(:,2)==b,1);
        if numel(parents_b)~=1||parents_b(1)~=a, continue, end
        if ~same_statement_group(groups(a),groups(b)), continue, end
        if groups(a).executions~=groups(b).executions, continue, end
        if ~elementwise_group(groups(a))||~elementwise_group(groups(b)), continue, end
        if any(feedback(:,1)==a|feedback(:,2)==a|feedback(:,1)==b|feedback(:,2)==b), continue, end
        groups(a)=combine_groups(groups(a),groups(b));
        [groups,edges,feedback]=remove_group(groups,edges,feedback,b,a);
        changed=true; break
    end
end
for k=1:numel(groups), groups(k).id=k; groups(k)=finalize_group(groups(k)); end
end

function tf=same_statement_group(a,b)
tf=strcmp(a.statement_key,b.statement_key)&&strcmp(a.call_path,b.call_path)&&strcmp(a.loop_id,b.loop_id);
end

function tf=elementwise_group(g)
allowed={'plus','minus','times','rdivide','ldivide','power','sqrt','abs'};
tf=~isempty(g.operators)&&all(ismember(g.operators,allowed));
end

function g=combine_groups(a,b)
g=a;
g.counts=add_counts(a.counts,b.counts);
g.operators=[a.operators b.operators];
g.output_variants=b.output_variants;
g.slot_start=min(a.slot_start,b.slot_start); g.slot_end=max(a.slot_end,b.slot_end);
g.first_raw_id=min(a.first_raw_id,b.first_raw_id);
g.last_raw_id=max(a.last_raw_id,b.last_raw_id);
g=finalize_group(g);
end

function [groups,edges,feedback]=remove_group(groups,edges,feedback,drop,keep)
% Redirect references to DROP into KEEP, delete DROP, and renumber IDs.
old_count=numel(groups);
edges(edges==drop)=keep; feedback(feedback==drop)=keep;
edges=unique(edges(edges(:,1)~=edges(:,2),:),'rows','stable');
feedback=unique(feedback,'rows','stable');
groups(drop)=[];
old_to_new=zeros(1,old_count);
next=0;
for old=1:old_count
    if old==drop, continue, end
    next=next+1; old_to_new(old)=next;
end
old_to_new(drop)=old_to_new(keep);
if ~isempty(edges)
    e1=arrayfun(@(x) old_to_new(x),edges(:,1));
    e2=arrayfun(@(x) old_to_new(x),edges(:,2));
    edges=[e1(:) e2(:)];
end
if ~isempty(feedback)
    f1=arrayfun(@(x) old_to_new(x),feedback(:,1));
    f2=arrayfun(@(x) old_to_new(x),feedback(:,2));
    feedback=[f1(:) f2(:)];
end
edges=unique(edges,'rows','stable'); feedback=unique(feedback,'rows','stable');
end

function info=empty_loop_info()
info=struct('id','','label','','line',0,'end_line',0);
end

function info=nearest_loop_context(frames,root_function)
info=empty_loop_info();
frames=frames_within_trace_root(frames,root_function);
if isempty(frames)||~isstruct(frames), return, end
% Use the nearest active loop inside the selected replay root. Frames above
% the selected top-level function belong to the original caller that merely
% invoked the measured region and must not create loop regions in its graph.
% Frames below the root are retained so a helper called from inside a loop
% in the selected function still inherits that loop context.
for k=1:numel(frames)
    candidate=active_loop_at(frames(k).file,frames(k).line);
    if ~isempty(candidate.id), info=candidate; return, end
end
end

function frames=frames_within_trace_root(frames,root_function)
% Keep only stack frames belonging to the selected function and its callees.
% dbstack is ordered from the arithmetic operation outward. The last frame
% whose leaf function name matches ROOT_FUNCTION is the replay boundary;
% frames above it are from the user's caller and are outside the analyzed
% function. Using the last match also behaves sensibly for recursive calls.
if isempty(frames)||~isstruct(frames), return, end
root_leaf=function_leaf_name(root_function);
names=cell(1,numel(frames));
for k=1:numel(frames), names{k}=function_leaf_name(frames(k).function_name); end
boundary=find(strcmpi(names,root_leaf),1,'last');
if isempty(boundary)
    % Defensive fallback: never import caller control flow when the root
    % frame cannot be identified. The innermost arithmetic frame is still
    % sufficient for direct operations in ordinary standalone functions.
    frames=frames(1);
else
    frames=frames(1:boundary);
end
end

function leaf=function_leaf_name(name)
leaf='';
if isempty(name), return, end
parts=regexp(char(name),'[./>]','split');
if ~isempty(parts), leaf=parts{end}; end
end

function info=active_loop_at(file,line_number)
info=empty_loop_info();
if isempty(file)||line_number<1||exist(file,'file')~=2, return, end
try, lines=regexp(fileread(file),'\r\n|\n|\r','split'); catch, return, end
if isempty(lines), return, end

function_start=enclosing_function_start(lines,line_number);
if function_start<1, function_start=1; end
stack=struct('type',{},'line',{},'label',{});
openers={'if','for','parfor','while','switch','try','spmd','arguments'};
last_line=min(line_number-1,numel(lines));
for k=function_start:last_line
    code=strip_comment(lines{k});
    statements=split_top_level_control_statements(code);
    for q=1:numel(statements)
        statement=strtrim(statements{q});
        token=regexp(statement,'^([A-Za-z]\w*)','tokens','once');
        if isempty(token), continue, end
        word=lower(token{1});
        if any(strcmp(word,openers))
            stack(end+1)=struct('type',word,'line',k,'label',statement); %#ok<AGROW>
        elseif strcmp(word,'end')&&~isempty(stack)
            stack(end)=[];
        end
    end
end
for k=numel(stack):-1:1
    if any(strcmp(stack(k).type,{'for','parfor','while'}))
        info.id=sprintf('%s|L%d',lower(file),stack(k).line);
        info.label=stack(k).label;
        info.line=stack(k).line;
        info.end_line=matching_control_end(lines,stack(k).line);
        return
    end
end
end

function line_number=enclosing_function_start(lines,target_line)
line_number=0;
for k=1:min(target_line,numel(lines))
    code=strtrim(strip_comment(lines{k}));
    if ~isempty(regexp(code,'^function(?:\s|$)','once'))
        line_number=k;
    end
end
end

function end_line=matching_control_end(lines,start_line)
end_line=0; depth=0;
openers={'if','for','parfor','while','switch','try','spmd','arguments'};
for k=start_line:numel(lines)
    code=strip_comment(lines{k});
    statements=split_top_level_control_statements(code);
    for q=1:numel(statements)
        statement=strtrim(statements{q});
        token=regexp(statement,'^([A-Za-z]\w*)','tokens','once');
        if isempty(token), continue, end
        word=lower(token{1});
        if any(strcmp(word,openers))
            depth=depth+1;
        elseif strcmp(word,'end')
            depth=depth-1;
            if depth==0, end_line=k; return, end
        end
    end
end
end

function pieces=split_top_level_control_statements(text)
% MATLAB permits compact blocks such as "for k=1:n, x=x+1; end". Split
% only top-level commas/semicolons, never delimiters inside (), [], {}, or
% strings, so those openers and their matching end are both observed.
pieces={}; start=1; dr=0; ds=0; dc=0; in_string=false; i=1;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''', i=i+2; continue, end
        prev=previous_nonspace(text(1:max(i-1,0)));
        if ~in_string&&~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}'))
            i=i+1; continue
        end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(', dr=dr+1;
        elseif ch==')', dr=max(0,dr-1);
        elseif ch=='[', ds=ds+1;
        elseif ch==']', ds=max(0,ds-1);
        elseif ch=='{', dc=dc+1;
        elseif ch=='}', dc=max(0,dc-1);
        elseif (ch==';'||ch==',')&&dr==0&&ds==0&&dc==0
            pieces{end+1}=text(start:i-1); %#ok<AGROW>
            start=i+1;
        end
    end
    i=i+1;
end
if start<=numel(text)
    pieces{end+1}=text(start:end); %#ok<AGROW>
elseif isempty(pieces)
    pieces={text};
end
end

function regions=empty_loop_regions()
regions=struct('id',{},'label',{},'line',{},'end_line',{}, ...
    'node_ids',{},'passes',{},'first_raw_id',{},'last_raw_id',{}, ...
    'file',{},'function_name',{},'call_path',{});
end

function regions=collect_loop_regions(groups)
regions=empty_loop_regions();
ids={groups.loop_id}; ids=unique(ids(~cellfun(@isempty,ids)),'stable');
for k=1:numel(ids)
    node_ids=find(strcmp({groups.loop_id},ids{k}));
    passes=max([groups(node_ids).executions]);
    regions(end+1)=struct('id',ids{k},'label',groups(node_ids(1)).loop_label, ...
        'line',groups(node_ids(1)).loop_line, ...
        'end_line',groups(node_ids(1)).loop_end_line, ...
        'node_ids',node_ids,'passes',passes, ...
        'first_raw_id',min([groups(node_ids).first_raw_id]), ...
        'last_raw_id',max([groups(node_ids).last_raw_id]), ...
        'file',groups(node_ids(1)).file, ...
        'function_name',groups(node_ids(1)).function_name, ...
        'call_path',groups(node_ids(1)).call_path); %#ok<AGROW>
end
end

function edges=add_loop_exit_barriers(groups,edges,regions)
% Preserve MATLAB control order without serializing ordinary independent
% statements. An operation lexically after a completed loop cannot execute
% before that loop finishes, even if it has no data dependency on the loop
% result. This prevents post-loop constants from being hoisted into the
% pre-loop parallel layer while retaining genuine parallel branches before
% the loop.
for r=1:numel(regions)
    ids=regions(r).node_ids;
    if isempty(ids)||regions(r).end_line<=0, continue, end

    internal_edges=edges(ismember(edges(:,1),ids)&ismember(edges(:,2),ids),:);
    sinks=ids;
    if ~isempty(internal_edges)
        sinks=setdiff(ids,unique(internal_edges(:,1)).','stable');
    end
    if isempty(sinks), continue, end

    candidates=[];
    for g=1:numel(groups)
        if ismember(g,ids), continue, end
        same_scope=strcmpi(groups(g).file,regions(r).file) && ...
            strcmp(groups(g).function_name,regions(r).function_name) && ...
            strcmp(groups(g).call_path,regions(r).call_path);
        after_loop=groups(g).first_raw_id>regions(r).last_raw_id && ...
            groups(g).line>=regions(r).end_line;
        if same_scope&&after_loop, candidates(end+1)=g; end %#ok<AGROW>
    end
    if isempty(candidates), continue, end

    % Connect only to entry nodes of the post-loop subgraph. Their
    % descendants inherit the barrier through normal data dependencies.
    candidate_parents=edges(ismember(edges(:,1),candidates)&ismember(edges(:,2),candidates),:);
    entries=candidates;
    if ~isempty(candidate_parents)
        entries=setdiff(candidates,unique(candidate_parents(:,2)).','stable');
    end
    for a=sinks(:).'
        for b=entries(:).'
            if a~=b, edges(end+1,:)=[a b]; end %#ok<AGROW>
        end
    end
end
if ~isempty(edges), edges=unique(edges,'rows','stable'); end
end

function groups=empty_group_nodes()
groups=repmat(empty_group_node(),1,0);
end

function g=empty_group_node()
g=struct('id',0,'first_raw_id',0,'last_raw_id',0, ...
    'file','','function_name','','line',0, ...
    'expression','','statement_key','','call_path','','operators',{{}}, ...
    'operator_text','','input_variants',{{}},'output_variants',{{}}, ...
    'dimension_text','','operation_text','','counts',zero_counts(), ...
    'count_text','','executions',0,'slot_start',1,'slot_end',1, ...
    'statement_slots',1,'loop_id','','loop_label','','loop_line',0, ...
    'loop_end_line',0,'depth',1,'x',0,'y',0);
end

function text=operator_sequence_text(ops)
if isempty(ops), text=''; return, end
pretty=cell(size(ops));
for k=1:numel(ops), pretty{k}=upper(ops{k}); end
text=strjoin(pretty,' -> ');
end

function text=group_dimension_text(g)
in=variant_text(g.input_variants); out=variant_text(g.output_variants);
if isempty(in), in='(input shape unavailable)'; end
if isempty(out), out='(output shape unavailable)'; end
text=sprintf('%s -> %s',in,out);
end

function text=variant_text(values)
if isempty(values), text=''; return, end
if numel(values)==1, text=values{1}; return, end
if numel(values)==2, text=sprintf('%s | %s',values{1},values{2}); return, end
text=sprintf('%s | %s | +%d variants',values{1},values{2},numel(values)-2);
end

function text=operation_text(c)
parts={};
if c.complex_mult>0, parts{end+1}='CMUL'; end %#ok<AGROW>
if c.complex_add>0, parts{end+1}='CADD'; end %#ok<AGROW>
if c.real_mult>0, parts{end+1}='RMUL'; end %#ok<AGROW>
if c.real_add>0, parts{end+1}='RADD'; end %#ok<AGROW>
if c.divisions>0, parts{end+1}='DIV'; end %#ok<AGROW>
if c.square_roots>0, parts{end+1}='SQRT'; end %#ok<AGROW>
if c.sort_comparisons>0 parts{end+1}='SORTCMP';end
text=strjoin(parts,' + ');
end

function text=count_text(c)
parts={};
if c.complex_mult>0, parts{end+1}=sprintf('%g CMUL',c.complex_mult); end %#ok<AGROW>
if c.complex_add>0, parts{end+1}=sprintf('%g CADD',c.complex_add); end %#ok<AGROW>
if c.real_mult>0, parts{end+1}=sprintf('%g RMUL',c.real_mult); end %#ok<AGROW>
if c.real_add>0, parts{end+1}=sprintf('%g RADD',c.real_add); end %#ok<AGROW>
if c.divisions>0, parts{end+1}=sprintf('%g DIV',c.divisions); end %#ok<AGROW>
if c.square_roots>0, parts{end+1}=sprintf('%g SQRT',c.square_roots); end %#ok<AGROW>
if c.sort_comparisons>0
    parts{end+1}=sprintf('%g SORTCMP',c.sort_comparisons);
end
text=strjoin(parts,' | ');
end

function [depth,max_width,critical]=graph_metrics(groups,edges)
if isempty(groups), depth=[]; max_width=0; critical=0; return, end
depth=ones(1,numel(groups)); weighted=max(1,[groups.executions]);
[~,order]=sort([groups.first_raw_id]);
for idx=order
    parents=edges(edges(:,2)==idx,1);
    if ~isempty(parents)
        depth(idx)=1+max(depth(parents));
        weighted(idx)=max(1,groups(idx).executions)+max(weighted(parents));
    end
end
max_width=0;
for d=1:max(depth), max_width=max(max_width,sum(depth==d)); end
critical=max(weighted);
end

function summaries=build_source_line_summary(groups)
summaries=struct('file',{},'function_name',{},'line',{},'expression',{}, ...
    'operation_text',{},'node_count',{},'max_executions',{},'counts_internal',{});
keys={};
for k=1:numel(groups)
    key=sprintf('%s|%s|%d|%s',lower(groups(k).file),lower(groups(k).function_name), ...
        groups(k).line,regexprep(strtrim(groups(k).expression),'\s+',' '));
    idx=find(strcmp(keys,key),1);
    if isempty(idx)
        keys{end+1}=key; %#ok<AGROW>
        idx=numel(keys);
        summaries(end+1)=struct('file',groups(k).file,'function_name',groups(k).function_name, ...
            'line',groups(k).line,'expression',groups(k).expression, ...
            'operation_text','','node_count',0,'max_executions',1, ...
            'counts_internal',zero_counts()); %#ok<AGROW>
        counts=zero_counts();
    else
        counts=summaries(idx).counts_internal;
    end
    counts=add_counts(counts,groups(k).counts);
    summaries(idx).counts_internal=counts;
    summaries(idx).node_count=summaries(idx).node_count+1;
    summaries(idx).max_executions=max(summaries(idx).max_executions,groups(k).executions);
end
for k=1:numel(summaries)
    summaries(k).operation_text=operation_text(summaries(k).counts_internal);
end
% Remove the temporary internal count field while preserving a homogeneous
% struct array.
if ~isempty(summaries), summaries=rmfield(summaries,'counts_internal'); end
[~,order]=sortrows([[summaries.line].' (1:numel(summaries)).']); summaries=summaries(order);
end

function print_report(stats)
fprintf('\nOPCOUNT REPORT\n');
fprintf('Function                : %s\n',stats.function_name);
fprintf('Source marker lines     : %d-%d\n',stats.source_region(1),stats.source_region(2));
fprintf('--------------------------------------------------\n');
fprintf('Complex multiplications : %15.0f\n',stats.complex_mult);
fprintf('Complex additions       : %15.0f\n',stats.complex_add);
fprintf('Real multiplications    : %15.0f\n',stats.real_mult);
fprintf('Real additions          : %15.0f\n',stats.real_add);
fprintf('Divisions               : %15.0f\n',stats.divisions);
fprintf('Square roots            : %15.0f\n',stats.square_roots);
fprintf('Sort comparisons        : %15.0f\n',stats.sort_comparisons);
fprintf('--------------------------------------------------\n');
fprintf('Real-equivalent mult.   : %15.0f\n',stats.real_equivalent_mult);
fprintf('Real-equivalent add.    : %15.0f\n',stats.real_equivalent_add);
fprintf('Measured original time  : %15.6f s\n',stats.elapsed_seconds);
fprintf('Dataflow operation nodes: %15.0f\n',stats.dataflow.operation_nodes);
fprintf('Maximum parallel width  : %15.0f\n',stats.dataflow.max_parallel_width);
fprintf('Critical-path stages    : %15.0f\n\n',stats.dataflow.critical_path_stages);

end

function plot_dataflow(flow,varargin)
if isempty(flow.nodes)
    error('opcount:EmptyGraph','The selected function executed no CMUL, CADD, RMUL, RADD, DIV, or SQRT operations.');
end
nodes=flow.nodes; edges=flow.edges; feedback=flow.feedback_edges;
[depth,max_width,~]=graph_metrics(nodes,edges);
[nodes,levels]=layout_nodes(nodes,edges,depth);
fig=figure('Name',['OPCOUNT Dataflow: ' flow.function_name], ...
    'NumberTitle','off','IntegerHandle','off','Tag','OPCOUNT_Dataflow', ...
    'Color','w','Position',[70 35 1580 max(900,115*levels)]);
ax=axes('Parent',fig,'Position',[0.025 0.055 0.700 0.885]); hold(ax,'on'); axis(ax,'off');
xlim(ax,[0 1]); ylim(ax,[0 1]);
loop_title=loop_title_text(flow.loop_regions);
title(ax,{sprintf('%s: coarse arithmetic dataflow',flow.function_name), ...
    sprintf('Top-to-bottom dependencies | %s | max width %d | weighted critical path %d stages', ...
    loop_title,max_width,flow.metrics.critical_path_stages)}, ...
    'Interpreter','none','FontWeight','bold','FontSize',12);

[node_w,node_h]=node_geometry(max_width,levels);
% Keep lexical pre-loop and post-loop operations visually outside the
% dashed loop boundary. This is a plotting-only adjustment: arithmetic
% counts, DAG dependencies, maximum width, and critical-path metrics are
% unchanged.
nodes=separate_loop_regions_for_display(nodes,edges,flow.loop_regions,node_w,node_h);
% Expand the plotting range when the display-only separation needs extra
% room above or below the nominal [0,1] layout.
y_values=[nodes.y];
ylim(ax,[min(0,min(y_values)-node_h/2-0.035) ...
    max(1,max(y_values)+node_h/2+0.035)]);
% Loop boundaries first.
for k=1:numel(flow.loop_regions)
    ids=flow.loop_regions(k).node_ids;
    ids=ids(ids<=numel(nodes)); if isempty(ids), continue, end
    x0=max(0.015,min([nodes(ids).x])-node_w/2-0.025);
    x1=min(0.985,max([nodes(ids).x])+node_w/2+0.025);
    y0=max(0.02,min([nodes(ids).y])-node_h/2-0.035);
    y1=min(0.97,max([nodes(ids).y])+node_h/2+0.035);
    rectangle(ax,'Position',[x0 y0 x1-x0 y1-y0],'Curvature',0.025, ...
        'LineStyle','--','LineWidth',1.15,'EdgeColor',[0.42 0.42 0.42]);
    label=sprintf('%s  (%d executions)',compact_expression(flow.loop_regions(k).label,55),flow.loop_regions(k).passes);
    text(ax,x0+0.008,y1-0.004,label,'Interpreter','none','FontWeight','bold', ...
        'FontSize',8.3,'Color',[0.28 0.28 0.28],'VerticalAlignment','top');
end
% Dependency arrows before nodes.
for k=1:size(edges,1), draw_graph_edge(ax,nodes(edges(k,1)),nodes(edges(k,2)),node_w,node_h); end
for k=1:size(feedback,1), draw_graph_feedback(ax,nodes(feedback(k,1)),nodes(feedback(k,2)),node_w); end
for k=1:numel(nodes), draw_graph_node(ax,nodes(k),node_w,node_h,flow.function_name,levels); end

% Right-side source panel in the requested order: source line first, then
% the operation classes on the next line.
panel_ax=axes('Parent',fig,'Position',[0.745 0.065 0.235 0.865]); axis(panel_ax,'off');
rectangle(panel_ax,'Position',[0 0 1 1],'Curvature',0.025, ...
    'FaceColor',[0.975 0.975 0.975],'EdgeColor',[0.75 0.75 0.75]);
text(panel_ax,0.05,0.965,'Counted source lines','FontWeight','bold','FontSize',11,'VerticalAlignment','top');
text(panel_ax,0.05,0.925,source_panel_text(flow.source_lines), ...
    'Interpreter','none','FontName','Courier New','FontSize',7.8, ...
    'VerticalAlignment','top','HorizontalAlignment','left');
text(panel_ax,0.05,0.315,'Operation color','FontWeight','bold','FontSize',10,'VerticalAlignment','top');
red_key=node_color(struct('complex_mult',1,'complex_add',0,'real_mult',0,'real_add',0,'divisions',0,'square_roots',0),1);
green_key=node_color(struct('complex_mult',0,'complex_add',1,'real_mult',0,'real_add',0,'divisions',0,'square_roots',0),1);
rectangle(panel_ax,'Position',[0.05 0.267 0.075 0.027],'Curvature',0.08, ...
    'FaceColor',red_key.face,'EdgeColor',red_key.edge,'LineWidth',1.1);
text(panel_ax,0.145,0.281,'MUL / DIV / SQRT','FontSize',8.2,'VerticalAlignment','middle');
rectangle(panel_ax,'Position',[0.05 0.226 0.075 0.027],'Curvature',0.08, ...
    'FaceColor',green_key.face,'EdgeColor',green_key.edge,'LineWidth',1.1);
text(panel_ax,0.145,0.240,'ADD / SUB only','FontSize',8.2,'VerticalAlignment','middle');
text(panel_ax,0.05,0.202,'Darker fill = more source executions','FontSize',8.0, ...
    'Color',[0.28 0.28 0.28],'VerticalAlignment','top');
text(panel_ax,0.05,0.160,'Counting policy','FontWeight','bold','FontSize',10,'VerticalAlignment','top');
policy={ ...
    'Counts follow the source literally;', ...
    'loop-invariant work is not hoisted.', ...
    '', ...
    'Guaranteed zero, diagonal, and', ...
    'triangular structure is exploited.', ...
    '', ...
    'Safe elementwise chains on one', ...
    'source statement are combined.'};
text(panel_ax,0.05,0.122,policy,'Interpreter','none','FontSize',8.1, ...
    'VerticalAlignment','top','Color',[0.28 0.28 0.28]);

if ~isempty(varargin)
    for k=1:2:numel(varargin)
        if k+1<=numel(varargin)&&strcmpi(varargin{k},'export')
            export_graph(fig,varargin{k+1});
        elseif k==1&&numel(varargin)==1&&(ischar(varargin{k})||isstring(varargin{k}))
            export_graph(fig,varargin{k});
        end
    end
end
end

function [nodes,levels]=layout_nodes(nodes,edges,depth)
levels=max(depth);
% Initial order is execution/source order. Two barycentric sweeps reduce
% crossings while retaining top-to-bottom dependencies.
orders=cell(1,levels);
for d=1:levels
    ids=find(depth==d); [~,q]=sort([nodes(ids).first_raw_id]); orders{d}=ids(q);
end
for sweep=1:2
    for d=2:levels
        ids=orders{d}; score=zeros(size(ids));
        prev=orders{d-1};
        for k=1:numel(ids)
            parents=edges(edges(:,2)==ids(k),1); pos=[];
            for p=parents(:).', j=find(prev==p,1); if ~isempty(j), pos(end+1)=j; end, end %#ok<AGROW>
            if isempty(pos), score(k)=numel(prev)+k; else, score(k)=mean(pos); end
        end
        [~,q]=sort(score); orders{d}=ids(q);
    end
    for d=levels-1:-1:1
        ids=orders{d}; score=zeros(size(ids)); nxt=orders{d+1};
        for k=1:numel(ids)
            children=edges(edges(:,1)==ids(k),2); pos=[];
            for c=children(:).', j=find(nxt==c,1); if ~isempty(j), pos(end+1)=j; end, end %#ok<AGROW>
            if isempty(pos), score(k)=numel(nxt)+k; else, score(k)=mean(pos); end
        end
        [~,q]=sort(score); orders{d}=ids(q);
    end
end
for d=1:levels
    ids=orders{d};
    if numel(ids)==1, xs=0.5; else, xs=linspace(0.10,0.90,numel(ids)); end
    y=0.93-(d-1)*(0.84/max(levels-1,1));
    for q=1:numel(ids), nodes(ids(q)).x=xs(q); nodes(ids(q)).y=y; nodes(ids(q)).depth=d; end
end
end

function nodes=separate_loop_regions_for_display(nodes,edges,regions,w,h)
% Move only visually overlapping non-loop nodes out of a loop's bounding
% box. The movement preserves their relative order and is based on lexical
% source position, so an initialization immediately before a loop is placed
% above it and all post-loop descendants are placed below it. The DAG itself
% is not modified.
if isempty(nodes)||isempty(regions), return, end
pad_x=0.025; pad_y=0.035; clearance=max(0.018,0.35*h);
for r=1:numel(regions)
    ids=regions(r).node_ids;
    ids=ids(ids>=1 & ids<=numel(nodes));
    if isempty(ids), continue, end

    x0=min([nodes(ids).x])-w/2-pad_x;
    x1=max([nodes(ids).x])+w/2+pad_x;
    y0=min([nodes(ids).y])-h/2-pad_y;
    y1=max([nodes(ids).y])+h/2+pad_y;

    outsiders=setdiff(1:numel(nodes),ids,'stable');
    intruders=[];
    for g=outsiders
        overlaps_x=(nodes(g).x+w/2)>x0 && (nodes(g).x-w/2)<x1;
        overlaps_y=(nodes(g).y+h/2)>y0 && (nodes(g).y-h/2)<y1;
        if overlaps_x&&overlaps_y, intruders(end+1)=g; end %#ok<AGROW>
    end
    if isempty(intruders), continue, end

    same_scope=false(1,numel(nodes));
    for g=1:numel(nodes)
        same_scope(g)=strcmpi(nodes(g).file,regions(r).file) && ...
            strcmp(nodes(g).function_name,regions(r).function_name) && ...
            strcmp(nodes(g).call_path,regions(r).call_path);
    end

    before=intruders(same_scope(intruders) & ...
        ([nodes(intruders).line] < regions(r).line));
    after=intruders(same_scope(intruders) & ...
        (regions(r).end_line>0) & ...
        ([nodes(intruders).line] >= regions(r).end_line));

    % Defensive fallback for generated/helper source where line ordering is
    % unavailable: use traced execution order relative to the loop body.
    unresolved=setdiff(intruders,[before after],'stable');
    for g=unresolved
        if nodes(g).last_raw_id<regions(r).first_raw_id
            before(end+1)=g; %#ok<AGROW>
        elseif nodes(g).first_raw_id>regions(r).last_raw_id
            after(end+1)=g; %#ok<AGROW>
        end
    end

    if ~isempty(before)
        % Only move the near-loop preamble; earlier independent branches
        % remain in their original top layers.
        pre=before([nodes(before).y]-h/2 < y1+clearance);
        if ~isempty(pre)
            current_bottom=min([nodes(pre).y])-h/2;
            desired_bottom=y1+clearance;
            shift=desired_bottom-current_bottom;
            for g=pre, nodes(g).y=nodes(g).y+shift; end
        end
    end

    if ~isempty(after)
        % Include post-loop descendants so their vertical dependency order
        % is preserved when the first post-loop layer is moved downward.
        post=after;
        frontier=after;
        while ~isempty(frontier)
            children=edges(ismember(edges(:,1),frontier),2).';
            children=children(~ismember(children,ids));
            children=children(same_scope(children));
            children=setdiff(unique(children,'stable'),post,'stable');
            post=[post children]; %#ok<AGROW>
            frontier=children;
        end
        current_top=max([nodes(post).y])+h/2;
        desired_top=y0-clearance;
        shift=desired_top-current_top;
        for g=post, nodes(g).y=nodes(g).y+shift; end
    end
end
end

function [w,h]=node_geometry(max_width,levels)
w=min(0.27,0.82/max(max_width,1)); h=min(0.092,0.68/max(levels,1));
h=max(h,0.052);
end

function draw_graph_node(ax,n,w,h,root_function,levels)
color=node_color(n.counts,n.executions);
rectangle(ax,'Position',[n.x-w/2 n.y-h/2 w h],'Curvature',0.06, ...
    'FaceColor',color.face,'EdgeColor',color.edge,'LineWidth',1.25);
expr=compact_expression(n.expression,66); if isempty(expr), expr='(source expression unavailable)'; end
if strcmp(short_function_name(n.function_name),short_function_name(root_function))
    tag=sprintf('[L%d]',n.line);
else
    tag=sprintf('[%s:L%d]',short_function_name(n.function_name),n.line);
end
if n.slot_start==n.slot_end
    step=sprintf('step %d/%d',n.slot_start,n.statement_slots);
else
    step=sprintf('steps %d-%d/%d',n.slot_start,n.slot_end,n.statement_slots);
end
opline=n.operation_text;
if n.statement_slots>1, opline=sprintf('%s   (%s)',opline,step); end
counts=n.count_text;
if n.executions>1, counts=sprintf('%s  (total for %d executions)',counts,n.executions); end
label={sprintf('%s  %s',expr,tag),opline,n.dimension_text,counts};
text(ax,n.x,n.y,label,'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','none','FontSize',max(6.6,8.3-levels/18),'Color',[0.08 0.08 0.08]);
end

function draw_graph_edge(ax,a,b,w,h)
dx=b.x-a.x; dy=b.y-a.y;
if abs(dy)>=abs(dx)
    if dy<0, p1=[a.x a.y-h/2]; p2=[b.x b.y+h/2];
    else, p1=[a.x a.y+h/2]; p2=[b.x b.y-h/2]; end
else
    if dx>0, p1=[a.x+w/2 a.y]; p2=[b.x-w/2 b.y];
    else, p1=[a.x-w/2 a.y]; p2=[b.x+w/2 b.y]; end
end
quiver(ax,p1(1),p1(2),p2(1)-p1(1),p2(2)-p1(2),0, ...
    'Color',[0.35 0.48 0.62],'LineWidth',1.0,'MaxHeadSize',0.18,'AutoScale','off');
end

function draw_graph_feedback(ax,a,b,w)
p1=[a.x+w/2 a.y]; xelbow=0.965; p3=[xelbow b.y]; p4=[b.x+w/2 b.y];
plot(ax,[p1(1) xelbow xelbow],[p1(2) p1(2) p3(2)],'--', ...
    'Color',[0.58 0.31 0.28],'LineWidth',1.0);
quiver(ax,p3(1),p3(2),p4(1)-p3(1),p4(2)-p3(2),0, ...
    'Color',[0.58 0.31 0.28],'LineStyle','--','LineWidth',1.0, ...
    'MaxHeadSize',0.20,'AutoScale','off');
text(ax,xelbow+0.006,(a.y+b.y)/2,'feedback','Rotation',90, ...
    'FontSize',7,'Color',[0.50 0.25 0.25],'HorizontalAlignment','center');
end

function text=source_panel_text(lines)
if isempty(lines), text={'Source information unavailable.'}; return, end
text={};
for k=1:numel(lines)
    expr=compact_expression(lines(k).expression,47);
    fun=short_function_name(lines(k).function_name);
    text{end+1}=sprintf('%s  [L%d, %s]',expr,lines(k).line,fun); %#ok<AGROW>
    text{end+1}=sprintf('  %s',lines(k).operation_text); %#ok<AGROW>
    extra=sprintf('  %d arithmetic node(s)',lines(k).node_count);
    if lines(k).max_executions>1, extra=sprintf('%s; up to %dx',extra,lines(k).max_executions); end
    text{end+1}=extra; text{end+1}=''; %#ok<AGROW>
end
end

function text=loop_title_text(regions)
if isempty(regions), text='no repeated loop body';
elseif numel(regions)==1, text=sprintf('%d loop passes',regions(1).passes);
else, text=sprintf('%d loop regions',numel(regions)); end
end

function export_graph(fig,filename)
if isstring(filename), filename=char(filename); end
if exist('exportgraphics','file')==2, exportgraphics(fig,filename,'Resolution',220);
else
    [~,~,ext]=fileparts(filename);
    if strcmpi(ext,'.pdf'), print(fig,filename,'-dpdf','-painters');
    else, print(fig,filename,'-dpng','-r220'); end
end
end

function c=node_color(counts,executions)
% Hardware-oriented color policy:
%   red family   -> multiplication, division, or square root is present
%   green family -> addition/subtraction only
% Darker saturation means the same source operation is executed repeatedly.
if nargin<2||isempty(executions), executions=1; end
repeat_level=min(1,log2(max(1,double(executions)))/3); % 1x,2x,4x,8x+
expensive=(counts.complex_mult>0)||(counts.real_mult>0)|| ...
    (counts.divisions>0)||(counts.square_roots>0);
add_only=~expensive&&((counts.complex_add>0)||(counts.real_add>0));
if expensive
    light_face=[1.00 0.94 0.94]; dark_face=[0.96 0.58 0.58];
    light_edge=[0.78 0.25 0.25]; dark_edge=[0.60 0.05 0.05];
elseif add_only
    light_face=[0.93 0.99 0.93]; dark_face=[0.56 0.86 0.60];
    light_edge=[0.20 0.58 0.28]; dark_edge=[0.04 0.38 0.12];
else
    light_face=[0.95 0.95 0.95]; dark_face=[0.78 0.78 0.78];
    light_edge=[0.45 0.45 0.45]; dark_edge=[0.25 0.25 0.25];
end
c.face=(1-repeat_level)*light_face+repeat_level*dark_face;
c.edge=(1-repeat_level)*light_edge+repeat_level*dark_edge;
end

function out=compact_expression(text,max_length)
if nargin<2, max_length=72; end
out=strtrim(regexprep(text,'\s+',' '));
if numel(out)>max_length, out=[out(1:max_length-3) '...']; end
end

function name=short_function_name(name)
parts=regexp(name,'[>.]','split'); name=parts{end};
end

% -------------------------------------------------------------------------
% Region parser: identifies exactly one top-level function call between the
% markers without requiring a wrapper around the user's function.
% -------------------------------------------------------------------------
function call = find_single_enclosed_call(file_name,start_line,stop_line)
source=fileread(file_name); lines=regexp(source,'\r\n|\n|\r','split');
first=start_line+1; last=stop_line-1;
if first>last || first>numel(lines), error('opcount:NoFunction','No function call exists between the markers.'); end
last=min(last,numel(lines));
statements=collect_statements(lines(first:last),first);
calls=struct('function_name',{},'input_expressions',{},'output_expressions',{},'statement',{},'line_number',{});
for k=1:numel(statements)
    c=parse_function_statement(statements(k).text,statements(k).line);
    if ~isempty(c), calls(end+1)=c; end %#ok<AGROW>
end
if isempty(calls), error('opcount:NoFunction','The region contains no top-level function invocation.'); end
if numel(calls)>1
    error('opcount:MultipleFunctions','The region contains %d top-level function calls (%s). Enclose one function at a time.',numel(calls),strjoin({calls.function_name},', '));
end
call=calls(1);
end

function statements=collect_statements(lines,first_line)
statements=struct('text',{},'line',{}); buffer=''; buffer_line=first_line;
dr=0; ds=0; dc=0; in_string=false;
for i=1:numel(lines)
    raw=strip_comment(lines{i}); trimmed=strtrim(raw);
    if isempty(trimmed)&&isempty(buffer), continue, end
    continuation=~isempty(regexp(trimmed,'\.\.\.\s*$','once'));
    if continuation, trimmed=regexprep(trimmed,'\.\.\.\s*$',''); end
    if isempty(buffer), buffer_line=first_line+i-1; else, buffer=[buffer ' ']; end %#ok<AGROW>
    buffer=[buffer trimmed]; %#ok<AGROW>
    [dr,ds,dc,in_string]=update_depths(trimmed,dr,ds,dc,in_string);
    if ~continuation&&dr==0&&ds==0&&dc==0&&~in_string
        pieces=split_top_level_semicolons(buffer);
        for j=1:numel(pieces)
            t=strtrim(pieces{j}); if ~isempty(t), statements(end+1)=struct('text',t,'line',buffer_line); end %#ok<AGROW>
        end
        buffer='';
    end
end
if ~isempty(strtrim(buffer)), statements(end+1)=struct('text',strtrim(buffer),'line',buffer_line); end
end

function line=strip_comment(line)
out=''; in_string=false; i=1;
while i<=numel(line)
    ch=line(i);
    if ch==''''
        if in_string&&i<numel(line)&&line(i+1)==''''
            out=[out ch line(i+1)]; i=i+2; continue %#ok<AGROW>
        end
        prev=previous_nonspace(out);
        if ~in_string&&~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}'))
            out=[out ch]; i=i+1; continue %#ok<AGROW>
        end
        in_string=~in_string; out=[out ch]; %#ok<AGROW>
    elseif ch=='%'&&~in_string, break
    else, out=[out ch]; %#ok<AGROW>
    end
    i=i+1;
end
line=out;
end

function ch=previous_nonspace(text)
idx=find(~isspace(text),1,'last'); if isempty(idx), ch=''; else, ch=text(idx); end
end

function [dr,ds,dc,in_string]=update_depths(text,dr,ds,dc,in_string)
i=1;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''', i=i+2; continue, end
        prev=previous_nonspace(text(1:max(i-1,0)));
        if ~in_string&&~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}')), i=i+1; continue, end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(',dr=dr+1;elseif ch==')',dr=dr-1;elseif ch=='[',ds=ds+1;elseif ch==']',ds=ds-1;elseif ch=='{',dc=dc+1;elseif ch=='}',dc=dc-1;end
    end
    i=i+1;
end
end

function pieces=split_top_level_semicolons(text)
pieces={}; start=1; dr=0;ds=0;dc=0;in_string=false;i=1;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''',i=i+2;continue,end
        prev=previous_nonspace(text(1:max(i-1,0)));
        if ~in_string&&~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}')),i=i+1;continue,end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(',dr=dr+1;elseif ch==')',dr=dr-1;elseif ch=='[',ds=ds+1;elseif ch==']',ds=ds-1;elseif ch=='{',dc=dc+1;elseif ch=='}',dc=dc-1;
        elseif ch==';'&&dr==0&&ds==0&&dc==0,pieces{end+1}=text(start:i-1);start=i+1;end %#ok<AGROW>
    end
    i=i+1;
end
if start<=numel(text),pieces{end+1}=text(start:end);elseif isempty(pieces),pieces={text};end
end

function call=parse_function_statement(statement,line_number)
call=[]; statement=strtrim(statement);
first_word=regexp(statement,'^[A-Za-z]\w*','match','once');
reserved={'if','elseif','else','end','for','parfor','while','switch','case','otherwise','try','catch','function','return','break','continue','global','persistent','classdef','properties','methods','events'};
if ~isempty(first_word)&&any(strcmp(first_word,reserved)),return,end
[eq,lhs]=top_level_assignment(statement);
if eq>0,rhs=strtrim(statement(eq+1:end));else,lhs='';rhs=statement;end
name_token=regexp(rhs,'^\s*([A-Za-z]\w*(?:\.[A-Za-z]\w*)*)\s*\(','tokens','once');
if isempty(name_token),return,end
open=find(rhs=='(',1,'first'); close=matching_parenthesis_local(rhs,open);
if close==0||~isempty(strtrim(rhs(close+1:end))),return,end
call=struct('function_name',name_token{1},'input_expressions',{split_top_level_commas(rhs(open+1:close-1))}, ...
    'output_expressions',{parse_lhs(lhs)},'statement',statement,'line_number',line_number);
end

function [pos,lhs]=top_level_assignment(text)
pos=0;lhs='';dr=0;ds=0;dc=0;in_string=false;i=1;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''',i=i+2;continue,end
        if ~in_string&&is_transpose_quote_local(text,i),i=i+1;continue,end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(',dr=dr+1;elseif ch==')',dr=dr-1;elseif ch=='[',ds=ds+1;elseif ch==']',ds=ds-1;elseif ch=='{',dc=dc+1;elseif ch=='}',dc=dc-1;
        elseif ch=='='&&dr==0&&ds==0&&dc==0
            prev=' ';nxt=' ';if i>1,prev=text(i-1);end;if i<numel(text),nxt=text(i+1);end
            if prev~='='&&prev~='<'&&prev~='>'&&prev~='~'&&nxt~='=',pos=i;lhs=strtrim(text(1:i-1));return,end
        end
    end
i=i+1;
end
end

function close=matching_parenthesis_local(text,open)
close=0;depth=0;in_string=false;i=open;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''',i=i+2;continue,end
        if ~in_string&&is_transpose_quote_local(text,i),i=i+1;continue,end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(',depth=depth+1;elseif ch==')',depth=depth-1;if depth==0,close=i;return,end,end
    end
i=i+1;
end
end

function tf=is_transpose_quote_local(text,pos)
prev='';for k=pos-1:-1:1,if ~isspace(text(k)),prev=text(k);break,end,end
tf=~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}'));
end

function outputs=parse_lhs(lhs)
if isempty(lhs),outputs={};return,end
lhs=strtrim(lhs);if startsWith(lhs,'[')&&endsWith(lhs,']'),outputs=split_top_level_commas(lhs(2:end-1));else,outputs={lhs};end
end

function parts=split_top_level_commas(text)
text=strtrim(text);if isempty(text),parts={};return,end
parts={};start=1;dr=0;ds=0;dc=0;in_string=false;i=1;
while i<=numel(text)
    ch=text(i);
    if ch==''''
        if in_string&&i<numel(text)&&text(i+1)=='''',i=i+2;continue,end
        prev=previous_nonspace(text(1:max(i-1,0)));
        if ~in_string&&~isempty(prev)&&(isstrprop(prev,'alphanum')||any(prev=='_])}')),i=i+1;continue,end
        in_string=~in_string;
    elseif ~in_string
        if ch=='(',dr=dr+1;elseif ch==')',dr=dr-1;elseif ch=='[',ds=ds+1;elseif ch==']',ds=ds-1;elseif ch=='{',dc=dc+1;elseif ch=='}',dc=dc-1;
        elseif ch==','&&dr==0&&ds==0&&dc==0,parts{end+1}=strtrim(text(start:i-1));start=i+1;end %#ok<AGROW>
    end
i=i+1;
end
parts{end+1}=strtrim(text(start:end));
end
%<OPCOUNT_TRACE_CLASS_BEGIN>
%|classdef opcount_trace_20260805j
%|    % Runtime-generated scalar tracing object used internally by opcount.m.
%|    % Each object carries a numeric value and the IDs of arithmetic nodes
%|    % that produced it. The class is generated automatically in tempdir.
%|
%|    properties
%|        value
%|        deps
%|        meta_only
%|        support
%|        unit_mask
%|    end
%|
%|    methods
%|        function obj = opcount_trace_20260805j(value,deps,meta_only,support,unit_mask)
%|            if nargin < 1, value = []; end
%|            if nargin < 2 || isempty(deps), deps = builtin('zeros',1,0); end
%|            if nargin < 3, meta_only = builtin('false'); end
%|            if nargin < 4 || isempty(support), support = builtin('true',size(value)); end
%|            if nargin < 5 || isempty(unit_mask), unit_mask = builtin('false',size(value)); end
%|            obj.value = value;
%|            obj.deps = unique(double(deps(:).'),'stable');
%|            obj.meta_only = logical(meta_only);
%|            % Normalize metadata masks at every construction boundary.  This
%|            % prevents a scalar or stale decomposition mask from reaching a
%|            % later elementwise operation with a different array shape.
%|            target_size = size(value);
%|            try
%|                support = opcount_trace_20260805j.expand_to(logical(support),target_size);
%|            catch
%|                support = builtin('true',target_size);
%|            end
%|            try
%|                unit_mask = opcount_trace_20260805j.expand_to(logical(unit_mask),target_size);
%|            catch
%|                unit_mask = builtin('false',target_size);
%|            end
%|            obj.support = logical(support);
%|            obj.unit_mask = logical(unit_mask) & obj.support;
%|        end
%|
%|        function out = plus(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av + bv;
%|            [sa,ua,sb,ub] = opcount_trace_20260805j.expanded_masks(a,b,value);
%|            overlap = sa & sb;
%|            counts = opcount_trace_20260805j.zero_counts();
%|            n = nnz(overlap);
%|            if ~isreal(av) && ~isreal(bv), counts.complex_add=n; else, counts.real_add=n; end
%|            support = sa | sb;
%|            unit = (sa & ~sb & ua) | (sb & ~sa & ub);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = minus(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av - bv;
%|            [sa,ua,sb,ub] = opcount_trace_20260805j.expanded_masks(a,b,value);
%|            overlap = sa & sb;
%|            counts = opcount_trace_20260805j.zero_counts();
%|            n = nnz(overlap);
%|            if ~isreal(av) && ~isreal(bv), counts.complex_add=n; else, counts.real_add=n; end
%|            support = sa | sb;
%|            unit = (sa & ~sb & ua) | (sb & ~sa & ub);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = uplus(a)
%|            out = a;
%|        end
%|
%|        function out = uminus(a)
%|            % Sign changes are intentionally not counted, matching the
%|            % coarse arithmetic policy used by OPCOUNT.
%|            out = opcount_trace_20260805j(-a.value,a.deps,a.meta_only,a.support,a.unit_mask);
%|        end
%|
%|        function out = times(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av .* bv;
%|            [sa,ua,sb,ub] = opcount_trace_20260805j.expanded_masks(a,b,value);
%|            active = sa & sb;
%|            mult_needed = active & ~(ua | ub);
%|            counts = opcount_trace_20260805j.mul_counts(nnz(mult_needed),0,av,bv);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,active,active & ua & ub);
%|        end
%|
%|        function out = mtimes(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av * bv;
%|            [products,sums,support,unit] = opcount_trace_20260805j.mtimes_structure(a,b,av,bv,value);
%|            counts = opcount_trace_20260805j.mul_counts(products,sums,av,bv);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = rdivide(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av ./ bv;
%|            [sa,ua,sb,ub] = opcount_trace_20260805j.expanded_masks(a,b,value); %#ok<ASGLU>
%|            counts = opcount_trace_20260805j.zero_counts();
%|            counts.divisions = nnz(sa & ~ub);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,sa,sa & ua & ub);
%|        end
%|
%|        function out = ldivide(a,b)
%|            out = rdivide(b,a);
%|        end
%|
%|        function out = mrdivide(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av / bv;
%|            if isscalar(bv)
%|                [sa,ua,sb,ub] = opcount_trace_20260805j.expanded_masks(a,b,value); %#ok<ASGLU>
%|                counts=opcount_trace_20260805j.zero_counts(); counts.divisions=nnz(sa & ~ub);
%|                support=sa; unit=sa & ua & ub;
%|            else
%|                [counts,support,unit]=opcount_trace_20260805j.solve_structure(b,a,bv,av,builtin('true'),value);
%|            end
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = mldivide(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av \ bv;
%|            if isscalar(av)
%|                [sb,ub,sa,ua] = opcount_trace_20260805j.expanded_masks(b,a,value); %#ok<ASGLU>
%|                counts=opcount_trace_20260805j.zero_counts(); counts.divisions=nnz(sb & ~ua);
%|                support=sb; unit=sb & ub & ua;
%|            else
%|                [counts,support,unit]=opcount_trace_20260805j.solve_structure(a,b,av,bv,builtin('false'),value);
%|            end
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = power(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av .^ bv;
%|            [counts,support,unit] = opcount_trace_20260805j.power_structure(a,av,bv,value,builtin('false'));
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = mpower(a,b)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            value = av ^ bv;
%|            [counts,support,unit] = opcount_trace_20260805j.power_structure(a,av,bv,value,builtin('true'));
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,support,unit);
%|        end
%|
%|        function out = sqrt(a)
%|            value = sqrt(a.value);
%|            counts = opcount_trace_20260805j.zero_counts();
%|            counts.square_roots = nnz(a.support & ~a.unit_mask);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts,a.support,a.unit_mask);
%|        end
%|
%|        function out = abs(a)
%|            value = abs(a.value);
%|            if isreal(a.value)
%|                out = opcount_trace_20260805j(value,a.deps,a.meta_only,a.support,a.unit_mask); return
%|            end
%|            active=a.support & ~a.unit_mask; n=nnz(active);
%|            counts = opcount_trace_20260805j.zero_counts();
%|            counts.real_mult = 2*n; counts.real_add = n; counts.square_roots = n;
%|            out = opcount_trace_20260805j.unary_result(a,value,counts,a.support,a.unit_mask);
%|        end
%|
%|        function out = norm(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = norm(a.value,args{:});
%|            active=nnz(a.support); mult_active=nnz(a.support & ~a.unit_mask);
%|            counts = opcount_trace_20260805j.zero_counts();
%|            if ~isreal(a.value)
%|                counts.real_mult = 2*mult_active; counts.real_add = max(2*active-1,0);
%|            else
%|                counts.real_mult = mult_active; counts.real_add = max(active-1,0);
%|            end
%|            counts.square_roots = double(active>0);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts,active>0,builtin('false'));
%|        end
%|
%|        function out = sum(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = sum(a.value,args{:});
%|            counts = opcount_trace_20260805j.zero_counts();
%|            nadd = max(numel(a.value)-numel(value),0);
%|            if isreal(a.value), counts.real_add = nadd; else, counts.complex_add = nadd; end
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function out = prod(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = prod(a.value,args{:});
%|            nmult = max(numel(a.value)-numel(value),0);
%|            counts = opcount_trace_20260805j.mul_counts(nmult,0,a.value,1);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function out = mean(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = mean(a.value,args{:});
%|            counts = opcount_trace_20260805j.zero_counts();
%|            nadd = max(numel(a.value)-numel(value),0);
%|            if isreal(a.value), counts.real_add = nadd; else, counts.complex_add = nadd; end
%|            counts.divisions = numel(value);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function out = dot(a,b,varargin)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = dot(av,bv,args{:});
%|            [sa,ua,sb,ub]=opcount_trace_20260805j.expanded_masks(a,b,av.*bv);
%|            active=sa&sb; products=nnz(active & ~(ua|ub));
%|            nterms=nnz(active); sums=max(nterms-numel(value),0);
%|            counts = opcount_trace_20260805j.mul_counts(products,sums,av,bv);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts,nterms>0,builtin('false'));
%|        end
%|
%|        function out = inv(a)
%|            value = inv(a.value);
%|            [counts,support,unit] = opcount_trace_20260805j.inv_structure(a);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts,support,unit);
%|        end
%|
%|        function out = pinv(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = pinv(a.value,args{:});
%|            [m,n] = size(a.value);
%|            q = min(m,n); r = max(m,n);
%|            mults = 2*r*q^2 + q^3;
%|            adds = mults;
%|            counts = opcount_trace_20260805j.zero_counts();
%|            if isreal(a.value)
%|                counts.real_mult = mults; counts.real_add = adds;
%|            else
%|                counts.complex_mult = mults; counts.complex_add = adds;
%|            end
%|            counts.divisions = q^2;
%|            counts.square_roots = q;
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function varargout = chol(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp = cell(1,max(nargout,1));
%|            [tmp{:}] = chol(a.value,args{:});
%|            n = size(a.value,1);
%|            counts = opcount_trace_20260805j.zero_counts();
%|            % Shahabuddin et al., Matrix Decomposition for Massive MIMO
%|            % Detection (NORCAS 2020), Algorithm 2 / Table III.  For a
%|            % complex n-by-n matrix: CM=(n^3-n)/6 and the real scaling of
%|            % off-diagonal complex values requires n(n-1) RMUL. Additions
%|            % below follow the literal recurrence (subtraction is ADD).
%|            if isreal(a.value)
%|                counts.real_mult = n*(n-1)*(n+4)/6;
%|                counts.real_add  = n*(n-1)*(n+1)/6;
%|            else
%|                counts.complex_mult = n*(n-1)*(n+1)/6;
%|                counts.real_mult    = n*(n-1);
%|                counts.complex_add  = n*(n-1)*(n-2)/6;
%|                counts.real_add     = n*(n-1)/2;
%|            end
%|            counts.divisions = n;       % reciprocals of diagonal elements
%|            counts.square_roots = n;
%|            id = opcount_trace_20260805j.record(counts,{a},tmp{1});
%|            lower_mode = builtin('false');
%|            for q=1:numel(args)
%|                if (ischar(args{q}) || (isstring(args{q}) && isscalar(args{q}))) && strcmpi(char(args{q}),'lower')
%|                    lower_mode = builtin('true');
%|                end
%|            end
%|            % Use guaranteed structural support from the requested
%|            % factorization, not incidental numerical zeros in this sample.
%|            % CHOL(A,'lower') is lower triangular; default CHOL is upper.
%|            if lower_mode
%|                support = tril(builtin('true',size(tmp{1})));
%|            else
%|                support = triu(builtin('true',size(tmp{1})));
%|            end
%|            tmp{1} = opcount_trace_20260805j(tmp{1},id,a.meta_only,support,builtin('false',size(tmp{1})));
%|            varargout = tmp(1:nargout);
%|        end
%|
%|        function varargout = qr(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp = cell(1,max(nargout,1));
%|            [tmp{:}] = qr(a.value,args{:});
%|            [m,n] = size(a.value); q = min(m,n);
%|            pairs = q*n - q*(q+1)/2;
%|            counts = opcount_trace_20260805j.zero_counts();
%|            % Classical Gram-Schmidt model from Shahabuddin et al. For a
%|            % square complex matrix this gives n^3 CMUL + 2n^2 RMUL, whose
%|            % real-equivalent multiplication count is n^2(4n+2).
%|            if isreal(a.value)
%|                counts.real_mult = 2*m*(q+pairs);
%|                counts.real_add  = q*(m-1) + pairs*(2*m-1);
%|            else
%|                counts.complex_mult = m*q + 2*m*pairs;
%|                counts.real_mult    = 2*m*q;
%|                counts.complex_add  = pairs*(2*m-1);
%|                counts.real_add     = q*(2*m-1);
%|            end
%|            counts.divisions = q; counts.square_roots = q;
%|            id = opcount_trace_20260805j.record(counts,{a},tmp{1});
%|            if nargout<=1
%|                % One-output QR returns R.
%|                support=triu(builtin('true',size(tmp{1}))) & (tmp{1}~=0 | logical(eye(size(tmp{1}))));
%|                tmp{1}=opcount_trace_20260805j(tmp{1},id,a.meta_only,support,builtin('false',size(tmp{1})));
%|            else
%|                % Q is generally dense/unitary; R is upper triangular.
%|                tmp{1}=opcount_trace_20260805j(tmp{1},id,a.meta_only,builtin('true',size(tmp{1})),builtin('false',size(tmp{1})));
%|                support=triu(builtin('true',size(tmp{2}))) & (tmp{2}~=0 | logical(eye(size(tmp{2}))));
%|                tmp{2}=opcount_trace_20260805j(tmp{2},id,a.meta_only,support,builtin('false',size(tmp{2})));
%|                if nargout>=3 && isnumeric(tmp{3}) && ismatrix(tmp{3}) && ~isvector(tmp{3})
%|                    sp=(tmp{3}~=0); un=sp & isreal(tmp{3}) & ((tmp{3}==1)|(tmp{3}==-1));
%|                    tmp{3}=opcount_trace_20260805j(tmp{3},id,a.meta_only,sp,un);
%|                end
%|            end
%|            varargout = tmp(1:nargout);
%|        end
%|
%|        function varargout = lu(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp = cell(1,max(nargout,1));
%|            [tmp{:}] = lu(a.value,args{:});
%|            [m,n]=size(a.value); q=min(m,n);
%|            products=0; divisions=0;
%|            for kk=1:q
%|                products=products+max(m-kk,0)*max(n-kk,0);
%|                divisions=divisions+max(m-kk,0);
%|            end
%|            counts=opcount_trace_20260805j.zero_counts();
%|            if isreal(a.value), counts.real_mult=products; counts.real_add=products;
%|            else, counts.complex_mult=products; counts.complex_add=products; end
%|            counts.divisions=divisions;
%|            id=opcount_trace_20260805j.record(counts,{a},tmp{1});
%|            vector_form=builtin('false');
%|            for qq=1:numel(args)
%|                if (ischar(args{qq}) || (isstring(args{qq})&&isscalar(args{qq}))) && strcmpi(char(args{qq}),'vector')
%|                    vector_form=builtin('true');
%|                end
%|            end
%|            if nargout<=1
%|                % Packed LU output: no standalone triangular guarantee.
%|                sp=(tmp{1}~=0); un=sp & isreal(tmp{1}) & ((tmp{1}==1)|(tmp{1}==-1));
%|                tmp{1}=opcount_trace_20260805j(tmp{1},id,a.meta_only,sp,un);
%|            else
%|                if nargout>=3
%|                    % With separate P, L is unit lower triangular.
%|                    sp=tril(builtin('true',size(tmp{1}))) & (tmp{1}~=0 | logical(eye(size(tmp{1}))));
%|                    un=builtin('false',size(tmp{1})); un(logical(eye(size(tmp{1}))))=builtin('true');
%|                else
%|                    % Two-output LU incorporates row permutation into L.
%|                    sp=(tmp{1}~=0); un=sp & isreal(tmp{1}) & ((tmp{1}==1)|(tmp{1}==-1));
%|                end
%|                tmp{1}=opcount_trace_20260805j(tmp{1},id,a.meta_only,sp,un);
%|                sp=triu(builtin('true',size(tmp{2}))) & (tmp{2}~=0 | logical(eye(size(tmp{2}))));
%|                tmp{2}=opcount_trace_20260805j(tmp{2},id,a.meta_only,sp,builtin('false',size(tmp{2})));
%|                for kk=3:min(nargout,numel(tmp))
%|                    if kk==5 && isnumeric(tmp{kk}) && ismatrix(tmp{kk})
%|                        sp=logical(eye(size(tmp{kk}))) & (tmp{kk}~=0); un=sp & isreal(tmp{kk}) & ((tmp{kk}==1)|(tmp{kk}==-1));
%|                        tmp{kk}=opcount_trace_20260805j(tmp{kk},id,a.meta_only,sp,un);
%|                    elseif ~vector_form && isnumeric(tmp{kk}) && ismatrix(tmp{kk}) && ~isvector(tmp{kk})
%|                        sp=(tmp{kk}~=0); un=sp & isreal(tmp{kk}) & ((tmp{kk}==1)|(tmp{kk}==-1));
%|                        tmp{kk}=opcount_trace_20260805j(tmp{kk},id,a.meta_only,sp,un);
%|                    end
%|                end
%|            end
%|            varargout=tmp(1:nargout);
%|        end
%|
%|        function varargout = ldl(a,varargin)
%|            args=opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp=cell(1,max(nargout,2));
%|            [tmp{:}]=ldl(a.value,args{:});
%|            n=size(a.value,1); counts=opcount_trace_20260805j.zero_counts();
%|            % Paper-matched complex multiplication model: real-equivalent
%|            % total (2n^3+12n^2-14n)/3 = 4*CMUL + RMUL.
%|            if isreal(a.value)
%|                counts.real_mult=n*(n-1)*(n+1)/6 + 2*n*(n-1);
%|                counts.real_add=n*(n-1)*(n+1)/6;
%|            else
%|                counts.complex_mult=n*(n-1)*(n+1)/6;
%|                counts.real_mult=4*n*(n-1);
%|                counts.complex_add=n*(n-1)*(n-2)/6;
%|                counts.real_add=n*(n-1)/2;
%|            end
%|            counts.divisions=n;
%|            id=opcount_trace_20260805j.record(counts,{a},tmp{1});
%|            vector_form=builtin('false');
%|            for qq=1:numel(args)
%|                if (ischar(args{qq}) || (isstring(args{qq})&&isscalar(args{qq}))) && strcmpi(char(args{qq}),'vector')
%|                    vector_form=builtin('true');
%|                end
%|            end
%|            if nargout>=1
%|                if nargout>=3
%|                    sp=tril(builtin('true',size(tmp{1}))) & (tmp{1}~=0 | logical(eye(size(tmp{1}))));
%|                    un=builtin('false',size(tmp{1})); un(logical(eye(size(tmp{1}))))=builtin('true');
%|                else
%|                    % Two-output LDL incorporates the symmetric permutation.
%|                    sp=(tmp{1}~=0); un=sp & isreal(tmp{1}) & ((tmp{1}==1)|(tmp{1}==-1));
%|                end
%|                tmp{1}=opcount_trace_20260805j(tmp{1},id,a.meta_only,sp,un);
%|            end
%|            if nargout>=2
%|                % D can contain 1x1 and 2x2 pivot blocks; retain its exact
%|                % block-diagonal support instead of falsely forcing diagonal.
%|                sp=(tmp{2}~=0); un=sp & isreal(tmp{2}) & ((tmp{2}==1)|(tmp{2}==-1));
%|                tmp{2}=opcount_trace_20260805j(tmp{2},id,a.meta_only,sp,un);
%|            end
%|            if nargout>=3 && ~vector_form && isnumeric(tmp{3}) && ismatrix(tmp{3}) && ~isvector(tmp{3})
%|                sp=(tmp{3}~=0); un=sp & isreal(tmp{3}) & ((tmp{3}==1)|(tmp{3}==-1));
%|                tmp{3}=opcount_trace_20260805j(tmp{3},id,a.meta_only,sp,un);
%|            end
%|            if nargout>=4 && isnumeric(tmp{4}) && ismatrix(tmp{4})
%|                sp=logical(eye(size(tmp{4}))) & (tmp{4}~=0); un=sp & isreal(tmp{4}) & ((tmp{4}==1)|(tmp{4}==-1));
%|                tmp{4}=opcount_trace_20260805j(tmp{4},id,a.meta_only,sp,un);
%|            end
%|            varargout=tmp(1:nargout);
%|        end
%|
%|        function out = fft(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = fft(a.value,args{:});
%|            n = opcount_trace_20260805j.transform_length(a.value,args);
%|            batches = max(numel(a.value)/max(n,1),1);
%|            stages = ceil(log2(max(n,1)));
%|            counts = opcount_trace_20260805j.zero_counts();
%|            counts.complex_mult = batches*(n/2)*stages;
%|            counts.complex_add = batches*n*stages;
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function out = ifft(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = ifft(a.value,args{:});
%|            n = opcount_trace_20260805j.transform_length(a.value,args);
%|            batches = max(numel(a.value)/max(n,1),1);
%|            stages = ceil(log2(max(n,1)));
%|            counts = opcount_trace_20260805j.zero_counts();
%|            counts.complex_mult = batches*(n/2)*stages;
%|            counts.complex_add = batches*n*stages;
%|            counts.divisions = numel(value);
%|            out = opcount_trace_20260805j.unary_result(a,value,counts);
%|        end
%|
%|        function out = conv(a,b,varargin)
%|            [av,bv] = opcount_trace_20260805j.values2(a,b);
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = conv(av,bv,args{:});
%|            nprod = numel(av)*numel(bv);
%|            nsum = max(nprod-numel(value),0);
%|            counts = opcount_trace_20260805j.mul_counts(nprod,nsum,av,bv);
%|            out = opcount_trace_20260805j.binary_result(a,b,value,counts);
%|        end
%|
%|        % Non-counted transforms: dependency and type are propagated.
%|        function out = conj(a), out = opcount_trace_20260805j(conj(a.value),a.deps,a.meta_only,a.support,a.unit_mask); end
%|        function out = real(a), out = opcount_trace_20260805j(real(a.value),a.deps,a.meta_only,a.support,a.unit_mask & isreal(a.value)); end
%|        function out = imag(a), if isreal(a.value), s=builtin('false',size(a.value)); else, s=a.support; end, out=opcount_trace_20260805j(imag(a.value),a.deps,a.meta_only,s,builtin('false',size(a.value))); end
%|        function out = transpose(a), out = opcount_trace_20260805j(a.value.',a.deps,a.meta_only,a.support.',a.unit_mask.'); end
%|        function out = ctranspose(a), out = opcount_trace_20260805j(a.value',a.deps,a.meta_only,a.support.',a.unit_mask.'); end
%|        function out = diag(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); v=diag(a.value,args{:}); s=diag(a.support,args{:}); u=diag(a.unit_mask,args{:}); out=opcount_trace_20260805j(v,a.deps,a.meta_only,s,u); end
%|        function out = triu(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(triu(a.value,args{:}),a.deps,a.meta_only,triu(a.support,args{:}),triu(a.unit_mask,args{:})); end
%|        function out = tril(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(tril(a.value,args{:}),a.deps,a.meta_only,tril(a.support,args{:}),tril(a.unit_mask,args{:})); end
%|        function out = reshape(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(reshape(a.value,args{:}),a.deps,a.meta_only,reshape(a.support,args{:}),reshape(a.unit_mask,args{:})); end
%|        function out = squeeze(a), out=opcount_trace_20260805j(squeeze(a.value),a.deps,a.meta_only,squeeze(a.support),squeeze(a.unit_mask)); end
%|        function out = permute(a,order), o=opcount_trace_20260805j.unwrap(order); out=opcount_trace_20260805j(permute(a.value,o),a.deps,a.meta_only,permute(a.support,o),permute(a.unit_mask,o)); end
%|        function out = ipermute(a,order), o=opcount_trace_20260805j.unwrap(order); out=opcount_trace_20260805j(ipermute(a.value,o),a.deps,a.meta_only,ipermute(a.support,o),ipermute(a.unit_mask,o)); end
%|        function out = repmat(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(repmat(a.value,args{:}),a.deps,a.meta_only,repmat(a.support,args{:}),repmat(a.unit_mask,args{:})); end
%|        function out = flip(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(flip(a.value,args{:}),a.deps,a.meta_only,flip(a.support,args{:}),flip(a.unit_mask,args{:})); end
%|        function out = fliplr(a), out=opcount_trace_20260805j(fliplr(a.value),a.deps,a.meta_only,fliplr(a.support),fliplr(a.unit_mask)); end
%|        function out = flipud(a), out=opcount_trace_20260805j(flipud(a.value),a.deps,a.meta_only,flipud(a.support),flipud(a.unit_mask)); end
%|        function out = circshift(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(circshift(a.value,args{:}),a.deps,a.meta_only,circshift(a.support,args{:}),circshift(a.unit_mask,args{:})); end
%|        function out = sign(a), out=opcount_trace_20260805j(sign(a.value),a.deps,a.meta_only,a.support,a.support); end
%|        function out = floor(a), out=opcount_trace_20260805j(floor(a.value),a.deps,a.meta_only); end
%|        function out = ceil(a), out=opcount_trace_20260805j(ceil(a.value),a.deps,a.meta_only); end
%|        function out = round(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(round(a.value,args{:}),a.deps,a.meta_only); end
%|        function out = fix(a), out=opcount_trace_20260805j(fix(a.value),a.deps,a.meta_only); end
%|        function out = mod(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=opcount_trace_20260805j(mod(av,bv),opcount_trace_20260805j.all_deps({a,b}),opcount_trace_20260805j.all_meta({a,b})); end
%|        function out = rem(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=opcount_trace_20260805j(rem(av,bv),opcount_trace_20260805j.all_deps({a,b}),opcount_trace_20260805j.all_meta({a,b})); end
%|        function out = exp(a), out=opcount_trace_20260805j(exp(a.value),a.deps,a.meta_only); end
%|        function out = log(a), out=opcount_trace_20260805j(log(a.value),a.deps,a.meta_only); end
%|        function out = log10(a), out=opcount_trace_20260805j(log10(a.value),a.deps,a.meta_only); end
%|        function out = sin(a), out=opcount_trace_20260805j(sin(a.value),a.deps,a.meta_only); end
%|        function out = cos(a), out=opcount_trace_20260805j(cos(a.value),a.deps,a.meta_only); end
%|        function out = tan(a), out=opcount_trace_20260805j(tan(a.value),a.deps,a.meta_only); end
%|        function out = asin(a), out=opcount_trace_20260805j(asin(a.value),a.deps,a.meta_only); end
%|        function out = acos(a), out=opcount_trace_20260805j(acos(a.value),a.deps,a.meta_only); end
%|        function out = atan(a), out=opcount_trace_20260805j(atan(a.value),a.deps,a.meta_only); end
%|        function out = atan2(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=opcount_trace_20260805j(atan2(av,bv),opcount_trace_20260805j.all_deps({a,b}),opcount_trace_20260805j.all_meta({a,b})); end
%|
%|        function varargout = min(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp = cell(1,max(nargout,1)); [tmp{:}] = min(a.value,args{:});
%|            tmp{1} = opcount_trace_20260805j(tmp{1},a.deps,a.meta_only);
%|            varargout = tmp(1:nargout);
%|        end
%|        function varargout = max(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            tmp = cell(1,max(nargout,1)); [tmp{:}] = max(a.value,args{:});
%|            tmp{1} = opcount_trace_20260805j(tmp{1},a.deps,a.meta_only);
%|            varargout = tmp(1:nargout);
%|        end
%|        function varargout = sort(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|
%|            % Perform actual MATLAB sorting
%|            tmp = cell(1,max(nargout,1));
%|            [tmp{:}] = sort(a.value,args{:});
%|
%|            % Number of elements being sorted
%|            n = numel(a.value);
%|
%|            % Estimated sorting comparison complexity
%|            counts = opcount_trace_20260805j.zero_counts();
%|            if n > 1
%|                counts.sort_comparisons = ceil(n*log2(n));
%|            end
%|
%|            % Record this SORT operation in OPCOUNT
%|            id = opcount_trace_20260805j.record( ...
%|                counts,{a},tmp{1});
%|
%|            if id == 0
%|                deps = a.deps;
%|            else
%|                deps = id;
%|            end
%|
%|            % Keep sorted output as traced object
%|            tmp{1} = opcount_trace_20260805j( ...
%|                tmp{1},deps,a.meta_only);
%|
%|            varargout = tmp(1:nargout);
%|        end
%|
%|        % Relational/logical operations determine control flow but are not
%|        % among the six requested arithmetic classes.
%|        function out = lt(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av<bv; end
%|        function out = le(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av<=bv; end
%|        function out = gt(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av>bv; end
%|        function out = ge(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av>=bv; end
%|        function out = eq(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av==bv; end
%|        function out = ne(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av~=bv; end
%|        function out = and(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av&bv; end
%|        function out = or(a,b), [av,bv]=opcount_trace_20260805j.values2(a,b); out=av|bv; end
%|        function out = not(a), out=~a.value; end
%|        function out = all(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=all(a.value,args{:}); end
%|        function out = any(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=any(a.value,args{:}); end
%|        function out = logical(a), out=logical(a.value); end
%|        function out = isreal(a), out=isreal(a.value); end
%|        function out = isnumeric(a), out=isnumeric(a.value); end
%|        function out = isfloat(a), out=isfloat(a.value); end
%|        function out = isinteger(a), out=isinteger(a.value); end
%|        function out = islogical(a), out=islogical(a.value); end
%|        function out = isscalar(a), out=isscalar(a.value); end
%|        function out = isvector(a), out=isvector(a.value); end
%|        function out = ismatrix(a), out=ismatrix(a.value); end
%|        function out = isrow(a), out=isrow(a.value); end
%|        function out = iscolumn(a), out=iscolumn(a.value); end
%|        function out = issparse(a), out=issparse(a.value); end
%|        function out = isempty(a), out=isempty(a.value); end
%|        function out = isfinite(a), out=isfinite(a.value); end
%|        function out = isnan(a), out=isnan(a.value); end
%|        function out = isinf(a), out=isinf(a.value); end
%|        function out = find(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=find(a.value,args{:}); end
%|
%|        function varargout = size(a,varargin)
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            if nargout <= 1
%|                varargout{1} = size(a.value,args{:});
%|            else
%|                tmp = cell(1,nargout); [tmp{:}] = size(a.value,args{:});
%|                varargout = tmp;
%|            end
%|        end
%|        function out = length(a), out=length(a.value); end
%|        function out = ndims(a), out=ndims(a.value); end
%|        function out = numel(a,varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=numel(a.value,args{:}); end
%|        function out = end(a,k,n), s=size(a.value); if k<n, out=s(k); else, out=prod(s(k:end)); end, end
%|        function out = double(a), out=double(a.value); end
%|        function out = single(a), out=single(a.value); end
%|        function out = char(a), out=char(a.value); end
%|
%|        function out = colon(a,b,varargin)
%|            if nargin==2
%|                out = colon(opcount_trace_20260805j.unwrap(a),opcount_trace_20260805j.unwrap(b));
%|            else
%|                out = colon(opcount_trace_20260805j.unwrap(a),opcount_trace_20260805j.unwrap(b),opcount_trace_20260805j.unwrap(varargin{1}));
%|            end
%|        end
%|
%|        function out = subsref(a,S)
%|            if strcmp(S(1).type,'.')
%|                out = builtin('subsref',a,S);
%|                return
%|            end
%|            value = builtin('subsref',a.value,S);
%|            out = opcount_trace_20260805j(value,a.deps,a.meta_only,builtin('subsref',a.support,S),builtin('subsref',a.unit_mask,S));
%|        end
%|
%|        function a = subsasgn(a,S,b)
%|            if strcmp(S(1).type,'.')
%|                a = builtin('subsasgn',a,S,b);
%|                return
%|            end
%|            bv = opcount_trace_20260805j.unwrap(b);
%|            [bs,bu]=opcount_trace_20260805j.masks(b);
%|            a.value = builtin('subsasgn',a.value,S,bv);
%|            a.support = builtin('subsasgn',a.support,S,bs);
%|            a.unit_mask = builtin('subsasgn',a.unit_mask,S,bu);
%|            if isa(b,'opcount_trace_20260805j')
%|                a.deps = unique([a.deps b.deps],'stable');
%|                a.meta_only = a.meta_only && b.meta_only;
%|            end
%|        end
%|
%|        function out = horzcat(varargin)
%|            out = opcount_trace_20260805j.cat_impl(2,varargin{:});
%|        end
%|        function out = vertcat(varargin)
%|            out = opcount_trace_20260805j.cat_impl(1,varargin{:});
%|        end
%|        function out = cat(dim,varargin)
%|            out = opcount_trace_20260805j.cat_impl(opcount_trace_20260805j.unwrap(dim),varargin{:});
%|        end
%|
%|        % Constructors dispatch here when at least one dimension/argument is
%|        % traced. They are non-counted and produce traced numeric arrays.
%|        function out = zeros(varargin), args=opcount_trace_20260805j.unwrap_args(varargin); v=builtin('zeros',args{:}); out=opcount_trace_20260805j(v,builtin('zeros',1,0),builtin('false'),builtin('false',size(v)),builtin('false',size(v))); end
%|        function out = ones(varargin), args=opcount_trace_20260805j.unwrap_args(varargin); v=builtin('ones',args{:}); out=opcount_trace_20260805j(v,builtin('zeros',1,0),builtin('false'),builtin('true',size(v)),builtin('true',size(v))); end
%|        function out = eye(varargin), args=opcount_trace_20260805j.unwrap_args(varargin); v=builtin('eye',args{:}); s=logical(v); out=opcount_trace_20260805j(v,builtin('zeros',1,0),builtin('false'),s,s); end
%|        function out = rand(varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(builtin('rand',args{:}),builtin('zeros',1,0),builtin('false')); end
%|        function out = randn(varargin), args=opcount_trace_20260805j.unwrap_args(varargin); out=opcount_trace_20260805j(builtin('randn',args{:}),builtin('zeros',1,0),builtin('false')); end
%|        function out = complex(a,b)
%|            if nargin<2, b=0; end
%|            [av,bv]=opcount_trace_20260805j.values2(a,b);
%|            value=complex(av,bv); [sa,ua,sb,ub]=opcount_trace_20260805j.expanded_masks(a,b,value); out=opcount_trace_20260805j(value,opcount_trace_20260805j.all_deps({a,b}),opcount_trace_20260805j.all_meta({a,b}),sa|sb,(sa&~sb&ua));
%|        end
%|    end
%|
%|    methods(Static)
%|        function out = trace_zeros(varargin)
%|            % Replay-only constructor used when a numeric ZEROS array later
%|            % receives traced indexed assignments.  It creates no arithmetic
%|            % node; it only preserves per-element support and dependencies.
%|            args = opcount_trace_20260805j.unwrap_args(varargin);
%|            value = builtin('zeros',args{:});
%|            mask = builtin('false',size(value));
%|            out = opcount_trace_20260805j(value,builtin('zeros',1,0),builtin('false'),mask,mask);
%|        end
%|    end
%|
%|    methods(Static, Access=private)
%|        function [av,bv] = values2(a,b)
%|            av = opcount_trace_20260805j.unwrap(a);
%|            bv = opcount_trace_20260805j.unwrap(b);
%|        end
%|
%|        function value = unwrap(x)
%|            if isa(x,'opcount_trace_20260805j'), value=x.value; else, value=x; end
%|        end
%|
%|        function args = unwrap_args(args)
%|            for k=1:numel(args), args{k}=opcount_trace_20260805j.unwrap(args{k}); end
%|        end
%|
%|        function deps = all_deps(args)
%|            deps=builtin('zeros',1,0);
%|            for k=1:numel(args)
%|                if isa(args{k},'opcount_trace_20260805j'), deps=[deps args{k}.deps]; end %#ok<AGROW>
%|            end
%|            deps=unique(deps,'stable');
%|        end
%|
%|        function tf = all_meta(args)
%|            seen=builtin('false'); tf=builtin('true');
%|            for k=1:numel(args)
%|                if isa(args{k},'opcount_trace_20260805j')
%|                    seen=builtin('true'); tf=tf && args{k}.meta_only;
%|                end
%|            end
%|            tf=seen && tf;
%|        end
%|
%|        function out = binary_result(a,b,value,counts,support,unit)
%|            if nargin<5, support=builtin('true',size(value)); end
%|            if nargin<6, unit=builtin('false',size(value)); end
%|            args={a,b}; meta=opcount_trace_20260805j.all_meta(args);
%|            if meta
%|                out=opcount_trace_20260805j(value,opcount_trace_20260805j.all_deps(args),builtin('true'),support,unit); return
%|            end
%|            id=opcount_trace_20260805j.record(counts,args,value,support,unit);
%|            if id==0, deps=opcount_trace_20260805j.all_deps(args); else, deps=id; end
%|            out=opcount_trace_20260805j(value,deps,builtin('false'),support,unit);
%|        end
%|
%|        function out = unary_result(a,value,counts,support,unit)
%|            if nargin<4, support=builtin('true',size(value)); end
%|            if nargin<5, unit=builtin('false',size(value)); end
%|            if a.meta_only
%|                out=opcount_trace_20260805j(value,a.deps,builtin('true'),support,unit); return
%|            end
%|            id=opcount_trace_20260805j.record(counts,{a},value,support,unit);
%|            if id==0, deps=a.deps; else, deps=id; end
%|            out=opcount_trace_20260805j(value,deps,builtin('false'),support,unit);
%|        end
%|
%|        function counts = zero_counts()
%|            counts=struct('complex_mult',0,'complex_add',0,'real_mult',0,'real_add',0,'divisions',0,'square_roots',0, 'sort_comparisons',0);
%|        end
%|
%|        function counts = mul_counts(products,sums,av,bv)
%|            counts=opcount_trace_20260805j.zero_counts();
%|            if ~isreal(av) && ~isreal(bv)
%|                counts.complex_mult=products; counts.complex_add=sums;
%|            elseif ~isreal(av) || ~isreal(bv)
%|                counts.real_mult=2*products; counts.real_add=2*sums;
%|            else
%|                counts.real_mult=products; counts.real_add=sums;
%|            end
%|        end
%|
%|        function [support,unit] = masks(x)
%|            if isa(x,'opcount_trace_20260805j')
%|                support=x.support; unit=x.unit_mask;
%|            else
%|                support=(x~=0);
%|                unit=isreal(x) & ((x==1)|(x==-1));
%|            end
%|            support=logical(support); unit=logical(unit)&support;
%|        end
%|
%|        function [sa,ua,sb,ub] = expanded_masks(a,b,value)
%|            [sa,ua]=opcount_trace_20260805j.masks(a); [sb,ub]=opcount_trace_20260805j.masks(b);
%|            target=builtin('true',size(value));
%|            try, sa=sa & target; ua=ua & target; catch, sa=opcount_trace_20260805j.expand_to(sa,size(value)); ua=opcount_trace_20260805j.expand_to(ua,size(value)); end
%|            try, sb=sb & target; ub=ub & target; catch, sb=opcount_trace_20260805j.expand_to(sb,size(value)); ub=opcount_trace_20260805j.expand_to(ub,size(value)); end
%|        end
%|
%|        function out = expand_to(mask,target_size)
%|            if isequal(size(mask),target_size), out=mask; return, end
%|            if isscalar(mask), out=repmat(mask,target_size); return, end
%|            sm=size(mask); nd=max(numel(sm),numel(target_size)); sm(end+1:nd)=1; ts=target_size; ts(end+1:nd)=1;
%|            reps=ones(1,nd);
%|            for q=1:nd
%|                if sm(q)==ts(q), reps(q)=1; elseif sm(q)==1, reps(q)=ts(q); else, error('opcount:MaskExpansion','Incompatible implicit expansion.'); end
%|            end
%|            out=repmat(mask,reps);
%|        end
%|
%|        function [products,sums,support,unit] = mtimes_structure(a,b,av,bv,value)
%|            [sa,ua]=opcount_trace_20260805j.masks(a); [sb,ub]=opcount_trace_20260805j.masks(b);
%|            if isscalar(av) || isscalar(bv)
%|                [sa,ua,sb,ub]=opcount_trace_20260805j.expanded_masks(a,b,value);
%|                active=sa&sb; products=nnz(active & ~(ua|ub)); sums=0;
%|                support=active; unit=active&ua&ub; return
%|            end
%|            m=size(av,1); k=size(av,2); n=size(bv,2);
%|            support=builtin('false',m,n); unit=builtin('false',m,n); products=0; sums=0;
%|            for ii=1:m
%|                for jj=1:n
%|                    active=sa(ii,1:k) & sb(1:k,jj).';
%|                    t=nnz(active); support(ii,jj)=t>0;
%|                    if t>0
%|                        free=ua(ii,1:k) | ub(1:k,jj).';
%|                        products=products+nnz(active & ~free);
%|                        sums=sums+max(t-1,0);
%|                        if t==1
%|                            qq=find(active,1); unit(ii,jj)=ua(ii,qq)&&ub(qq,jj);
%|                        end
%|                    end
%|                end
%|            end
%|        end
%|
%|        function kind = structure_kind(support,unit)
%|            support=logical(support); unit=logical(unit);
%|            % Defensive normalization: structural and unit masks must always
%|            % describe the same array before any indexed classification.
%|            if ~isequal(size(unit),size(support))
%|                try, unit=opcount_trace_20260805j.expand_to(unit,size(support));
%|                catch, unit=builtin('false',size(support)); end
%|            end
%|            unit=unit & support;
%|            if isempty(support)||~any(support(:)), kind='zero'; return, end
%|            if isscalar(support), if unit, kind='unit'; else, kind='scalar'; end, return, end
%|            [m,n]=size(support); if m~=n, kind='rectangular'; return, end
%|            I=logical(eye(n));
%|            if all(~support(~I)) && all(unit(I)), kind='identity'; return, end
%|            if all(~support(~I)), kind='diagonal'; return, end
%|            if all(~support(~triu(true(n),1))), kind='strict_upper'; return, end
%|            if all(~support(~triu(true(n),0))), kind='upper_triangular'; return, end
%|            if all(~support(~tril(true(n),-1))), kind='strict_lower'; return, end
%|            if all(~support(~tril(true(n),0))), kind='lower_triangular'; return, end
%|            % LDL can return 1x1/2x2 pivot blocks. Recognize a nonoverlapping
%|            % block diagonal rather than falling back to a dense matrix.
%|            band=I;
%|            if n>1
%|                first_offdiag=diag(builtin('true',n-1,1),1);
%|                band=band | first_offdiag | first_offdiag.';
%|            end
%|            if all(~support(~band))
%|                valid=builtin('true'); ii=1;
%|                while ii<=n
%|                    if ii<n && (support(ii,ii+1)||support(ii+1,ii))
%|                        if ii+1<n && (support(ii+1,ii+2)||support(ii+2,ii+1)), valid=builtin('false'); break, end
%|                        ii=ii+2;
%|                    else
%|                        ii=ii+1;
%|                    end
%|                end
%|                if valid, kind='block_diagonal_2x2'; return, end
%|            end
%|            % Two-output LU/LDL incorporates a row permutation into L. A
%|            % row-permuted lower factor has prefix support with row counts
%|            % 1:n in a permuted order; permutation itself is wiring/indexing.
%|            rc=sum(support,2); lower_rows=builtin('true'); upper_rows=builtin('true');
%|            for ii=1:n
%|                c=rc(ii);
%|                if c<1 || c>n || any(support(ii,c+1:end)), lower_rows=builtin('false'); end
%|                if c<1 || c>n || any(support(ii,1:n-c)), upper_rows=builtin('false'); end
%|            end
%|            if isequal(sort(rc(:)).',(1:n))
%|                if lower_rows, kind='row_permuted_lower'; return, end
%|                if upper_rows, kind='row_permuted_upper'; return, end
%|            end
%|            if all(~support(I)), kind='off_diagonal'; return, end
%|            kind='dense';
%|        end
%|
%|        function [counts,support,unit] = inv_structure(a)
%|            n=size(a.value,1); counts=opcount_trace_20260805j.zero_counts();
%|            kind=opcount_trace_20260805j.structure_kind(a.support,a.unit_mask);
%|            support=a.support; unit=a.unit_mask;
%|            if strcmp(kind,'identity'), return, end
%|            if strcmp(kind,'diagonal')
%|                counts.divisions=nnz(diag(a.support) & ~diag(a.unit_mask)); return
%|            end
%|            if any(strcmp(kind,{'upper_triangular','lower_triangular'}))
%|                % Direct count for the supplied triangular-inverse recurrence:
%|                %   A(j,j) = 1/U(j,j)
%|                %   t = sum_{k=i}^{j-1} A(i,k)U(k,j)
%|                %   A(i,j) = -t*A(j,j)
%|                % Guaranteed zeros are skipped. The first active product seeds
%|                % the accumulator, so a d-term inner product uses d products
%|                % and d-1 additions, followed by one reciprocal scaling product.
%|                if contains(kind,'upper'), support=triu(true(n)); else, support=tril(true(n)); end
%|                unit=false(n);
%|                dot_mults=n*(n-1)*(n+1)/6;
%|                scale_mults=n*(n-1)/2;
%|                mults=dot_mults+scale_mults;
%|                adds=n*(n-1)*(n-2)/6;
%|                if isreal(a.value), counts.real_mult=mults; counts.real_add=adds; else, counts.complex_mult=mults; counts.complex_add=adds; end
%|                counts.divisions=n; return
%|            end
%|            if any(strcmp(kind,{'strict_upper','strict_lower'}))
%|                % A strictly triangular square matrix is singular. MATLAB INV
%|                % will fail before this model is used, but keep a conservative
%|                % dense fallback here for completeness.
%|                support=true(size(a.value)); unit=false(size(a.value));
%|            end
%|            support=true(size(a.value)); unit=false(size(a.value));
%|            if isreal(a.value), counts.real_mult=n^3; counts.real_add=n^3; else, counts.complex_mult=n^3; counts.complex_add=n^3; end
%|            counts.divisions=n^2;
%|        end
%|
%|        function [counts,support,unit] = solve_structure(Aobj,Bobj,A,B,right_division,value)
%|            if right_division
%|                % X=B/A is evaluated as a transposed left solve for cost/structure.
%|                At=A.'; Bt=B.'; As=Aobj.support.'; Au=Aobj.unit_mask.';
%|                [Bs,Bu]=opcount_trace_20260805j.masks(Bobj);
%|                temp=opcount_trace_20260805j(At,[],true,As,Au);
%|                tempB=opcount_trace_20260805j(Bt,[],true,Bs.',Bu.');
%|                [counts,st,ut]=opcount_trace_20260805j.solve_structure(temp,tempB,At,Bt,false,value.');
%|                support=st.'; unit=ut.'; return
%|            end
%|            n=size(A,1); p=size(B,2); counts=opcount_trace_20260805j.zero_counts();
%|            [bs,bu]=opcount_trace_20260805j.masks(Bobj); kind=opcount_trace_20260805j.structure_kind(Aobj.support,Aobj.unit_mask);
%|            if strcmp(kind,'identity'), support=bs; unit=bu; return, end
%|            if strcmp(kind,'diagonal')
%|                support=bs; unit=false(size(bs));
%|                nonunit_diag=~diag(Aobj.unit_mask);
%|                counts.divisions=sum(sum(bs & repmat(nonunit_diag,1,p))); return
%|            end
%|            if any(strcmp(kind,{'upper_triangular','lower_triangular','row_permuted_lower','row_permuted_upper'}))
%|                mults=p*n*(n-1)/2; adds=mults; support=true(size(value)); unit=false(size(value));
%|                if isreal(A)&&isreal(B), counts.real_mult=mults; counts.real_add=adds; else, counts.complex_mult=mults; counts.complex_add=adds; end
%|                counts.divisions=p*n; return
%|            end
%|            if strcmp(kind,'block_diagonal_2x2')
%|                support=bs; unit=false(size(bs)); mults=0; adds=0; divs=0; ii=1;
%|                while ii<=n
%|                    if ii<n && (Aobj.support(ii,ii+1)||Aobj.support(ii+1,ii))
%|                        % One 2x2 elimination is shared across RHS columns;
%|                        % each RHS then needs forward/back substitution.
%|                        mults=mults+1+2*p; adds=adds+1+2*p; divs=divs+1+2*p; ii=ii+2;
%|                    else
%|                        divs=divs+p; ii=ii+1;
%|                    end
%|                end
%|                if isreal(A)&&isreal(B), counts.real_mult=mults; counts.real_add=adds; else, counts.complex_mult=mults; counts.complex_add=adds; end
%|                counts.divisions=divs; return
%|            end
%|            mults=n^3+n^2*max(p,1); adds=mults; support=true(size(value)); unit=false(size(value));
%|            if ~isreal(A)||~isreal(B), counts.complex_mult=mults; counts.complex_add=adds; else, counts.real_mult=mults; counts.real_add=adds; end
%|            counts.divisions=n^2+n*max(p,1);
%|        end
%|
%|        function [counts,support,unit] = power_structure(a,base,exponent,value,is_matrix)
%|            counts=opcount_trace_20260805j.zero_counts();
%|            if ~isscalar(exponent)||~isreal(exponent)||exponent<0||exponent~=floor(exponent)
%|                support=true(size(value)); unit=false(size(value)); return
%|            end
%|            p=double(exponent);
%|            if p==0
%|                if is_matrix && ~isscalar(base), support=logical(eye(size(base,1))); unit=support; else, support=true(size(value)); unit=true(size(value)); end
%|                return
%|            elseif p==1
%|                support=a.support; unit=a.unit_mask; return
%|            end
%|            if ~is_matrix || isscalar(base)
%|                support=a.support; unit=a.unit_mask;
%|                products=(p-1)*nnz(a.support & ~a.unit_mask);
%|                if isreal(base), counts.real_mult=products; else, counts.complex_mult=products; end
%|                return
%|            end
%|            support=a.support; unit=a.unit_mask;
%|            current=opcount_trace_20260805j(base,[],true,support,unit);
%|            for q=2:p
%|                [prodq,sumq,support,unit]=opcount_trace_20260805j.mtimes_structure(current,a,current.value,base,current.value*base);
%|                cq=opcount_trace_20260805j.mul_counts(prodq,sumq,current.value,base); counts=opcount_trace_20260805j.add_count_structs(counts,cq);
%|                current=opcount_trace_20260805j(current.value*base,[],true,support,unit);
%|            end
%|        end
%|
%|        function c = add_count_structs(a,b)
%|            c=a; f=fieldnames(c); for q=1:numel(f), c.(f{q})=c.(f{q})+b.(f{q}); end
%|        end
%|
%|        function n = transform_length(value,args)
%|            if ~isempty(args) && isnumeric(args{1}) && isscalar(args{1})
%|                n=args{1}; return
%|            end
%|            s=size(value); d=find(s~=1,1,'first'); if isempty(d), d=1; end
%|            n=s(d);
%|        end
%|
%|        function out = cat_impl(dim,varargin)
%|            vals=cell(size(varargin)); ss=cell(size(varargin)); uu=cell(size(varargin));
%|            for k=1:numel(varargin), vals{k}=opcount_trace_20260805j.unwrap(varargin{k}); [ss{k},uu{k}]=opcount_trace_20260805j.masks(varargin{k}); end
%|            value=cat(dim,vals{:});
%|            out=opcount_trace_20260805j(value,opcount_trace_20260805j.all_deps(varargin),opcount_trace_20260805j.all_meta(varargin),cat(dim,ss{:}),cat(dim,uu{:}));
%|        end
%|
%|        function id = record(counts,args,value,support,unit)
%|            global OPCOUNT_RUNTIME_V28
%|            if nargin<4, support=builtin('true',size(value)); end
%|            if nargin<5, unit=builtin('false',size(value)); end
%|            id=0;
%|            if isempty(OPCOUNT_RUNTIME_V28) || ~isstruct(OPCOUNT_RUNTIME_V28) || ...
%|                    ~isfield(OPCOUNT_RUNTIME_V28,'active') || ~OPCOUNT_RUNTIME_V28.active
%|                return
%|            end
%|            names=fieldnames(counts); total=0;
%|            for k=1:numel(names), total=total+counts.(names{k}); end
%|            if total==0, return, end
%|            context=opcount_trace_20260805j.context();
%|            node=struct();
%|            node.id=OPCOUNT_RUNTIME_V28.next_id;
%|            node.parents=opcount_trace_20260805j.all_deps(args);
%|            node.counts=counts;
%|            node.file=context.file;
%|            node.function_name=context.function_name;
%|            node.line=context.line;
%|            node.expression=context.expression;
%|            node.operator=context.operator;
%|            node.stack_frames=context.frames;
%|            node.input_dimensions=cell(1,numel(args));
%|            for k=1:numel(args), [ss,uu]=opcount_trace_20260805j.masks(args{k}); node.input_dimensions{k}=opcount_trace_20260805j.describe(opcount_trace_20260805j.unwrap(args{k}),ss,uu); end
%|            node.output_dimensions=opcount_trace_20260805j.describe(value,support,unit);
%|            OPCOUNT_RUNTIME_V28.nodes(end+1)=node;
%|            OPCOUNT_RUNTIME_V28.next_id=OPCOUNT_RUNTIME_V28.next_id+1;
%|            id=node.id;
%|        end
%|
%|        function text = describe(value,support,unit)
%|            if nargin<2, support=builtin('true',size(value)); end
%|            if nargin<3, unit=builtin('false',size(value)); end
%|            s=size(value); parts=arrayfun(@num2str,s,'UniformOutput',builtin('false'));
%|            if isreal(value), kind='real'; else, kind='complex'; end
%|            structure=opcount_trace_20260805j.structure_kind(support,unit);
%|            if any(strcmp(structure,{'dense','rectangular','scalar'})), text=sprintf('[%s] %s',strjoin(parts,'x'),kind);
%|            else, text=sprintf('[%s] %s %s',strjoin(parts,'x'),kind,strrep(structure,'_',' ')); end
%|        end
%|
%|        function context = context()
%|            context=struct('file','','function_name','','line',0,'expression','', ...
%|                'operator','','frames',struct('file',{},'function_name',{},'line',{}));
%|            stack=dbstack('-completenames');
%|            helper_names={'context','record','binary_result','unary_result'};
%|            for k=2:numel(stack)
%|                f=stack(k).file;
%|                if isempty(f), continue, end
%|                [~,name]=fileparts(f);
%|                if strcmpi(name,'opcount_trace_20260805j')
%|                    if isempty(context.operator)
%|                        candidate=stack(k).name;
%|                        pieces=regexp(candidate,'[./>]','split'); candidate=lower(pieces{end});
%|                        if ~any(strcmp(candidate,helper_names)), context.operator=candidate; end
%|                    end
%|                    continue
%|                elseif strcmpi(name,'opcount')
%|                    continue
%|                end
%|                [mapped_file,mapped_line,mapped_name]=opcount_trace_20260805j.map_source_location(f,stack(k).line,stack(k).name);
%|                frame=struct('file',mapped_file,'function_name',mapped_name,'line',mapped_line);
%|                context.frames(end+1)=frame;
%|                if isempty(context.file)
%|                    context.file=mapped_file; context.function_name=mapped_name; context.line=mapped_line;
%|                    context.expression=opcount_trace_20260805j.source_statement(mapped_file,mapped_line);
%|                end
%|            end
%|        end
%|
%|        function [file_out,line_out,name_out] = map_source_location(file_in,line_in,name_in)
%|            global OPCOUNT_RUNTIME_V28
%|            file_out=file_in; line_out=line_in; name_out=name_in;
%|            try
%|                if isfield(OPCOUNT_RUNTIME_V28,'source_maps') && ~isempty(OPCOUNT_RUNTIME_V28.source_maps)
%|                    maps=OPCOUNT_RUNTIME_V28.source_maps;
%|                    for q=1:numel(maps)
%|                        if strcmpi(file_in,maps(q).generated_file) && line_in>=maps(q).generated_start
%|                            file_out=maps(q).original_file;
%|                            line_out=maps(q).original_start + (line_in-maps(q).generated_start);
%|                            parts=regexp(name_in,'[>/]','split');
%|                            name_out=parts{end};
%|                            return
%|                        end
%|                    end
%|                end
%|            catch
%|                file_out=file_in; line_out=line_in; name_out=name_in;
%|            end
%|        end
%|
%|        function text = source_statement(file,line_number)
%|            global OPCOUNT_RUNTIME_V28
%|            text='';
%|            try
%|                if ~isfield(OPCOUNT_RUNTIME_V28,'source_cache') || isempty(OPCOUNT_RUNTIME_V28.source_cache)
%|                    OPCOUNT_RUNTIME_V28.source_cache=containers.Map('KeyType','char','ValueType','any');
%|                end
%|                if isKey(OPCOUNT_RUNTIME_V28.source_cache,file)
%|                    lines=OPCOUNT_RUNTIME_V28.source_cache(file);
%|                else
%|                    lines=regexp(fileread(file),'\r\n|\n|\r','split');
%|                    OPCOUNT_RUNTIME_V28.source_cache(file)=lines;
%|                end
%|                if line_number<1 || line_number>numel(lines), return, end
%|                first=line_number;
%|                while first>1 && ~isempty(regexp(strtrim(lines{first-1}),'\.\.\.\s*$','once')), first=first-1; end
%|                last=line_number;
%|                while last<numel(lines) && ~isempty(regexp(strtrim(lines{last}),'\.\.\.\s*$','once')), last=last+1; end
%|                pieces=lines(first:last);
%|                for q=1:numel(pieces)
%|                    pieces{q}=regexprep(pieces{q},'%.*$','');
%|                    pieces{q}=regexprep(pieces{q},'\.\.\.\s*$','');
%|                    pieces{q}=strtrim(pieces{q});
%|                end
%|                text=strtrim(strjoin(pieces,' '));
%|            catch
%|                text='';
%|            end
%|        end
%|    end
%|end
%<OPCOUNT_TRACE_CLASS_END>
