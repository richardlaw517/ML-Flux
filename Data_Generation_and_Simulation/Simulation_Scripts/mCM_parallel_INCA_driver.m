% most up to date for parfor (chunked save/append + robust Constant)
function mCM_parallel_INCA_driver(numSetsOverride)
% mCM_parallel_driver (process-pool, worker-local model, no I/O in loop)
%
% Usage:
%   CCM_Nature_parallel_driver            % runs all Flux sets
%   CCM_Nature_parallel_driver(200)       % run first 200 flux sets
%
% What this script does:
%   • Builds the reaction network + MS data (specs below)
%   • Loads Flux sets from fluxfile.mat (expects variable "fluxset")
%   • Creates per-worker model copies via parallel.pool.Constant
%   • Runs simulate(m) for each flux set in parallel (parfor), in CHUNKS
%   • Overwrites sim_raw_chunk.mat per chunk; appends human-readable vals to raw_vals_named1.txt
%
% Notes for robust parallel execution:
%   • Requires a PROCESS-BASED pool (not threads). We'll ensure that.
%   • All helpers live in this function file as SUBFUNCTIONS so workers see them.
%   • No file I/O inside parfor. We collect in memory and save once per chunk.
%   • Each worker gets an independent model (no shared handle state).
% -------------------------------------------------------------------------

if nargin < 1, numSetsOverride = []; end

%% ---------- Model: reaction network ----------
reactionSpec = {
'CO20 (A) -> CO2 (A)';
'GLC12 (ABCDEF) -> GLC (ABCDEF)';
'GLC (ABCDEF) -> G6P (ABCDEF)';
'G6P (ABCDEF) <-> F6P (ABCDEF)';
'F6P (ABCDEF) -> FBP (ABCDEF)';
'FBP (ABCDEF) <-> DHAP (CBA) + GAP (DEF)';
'DHAP (ABC) <-> GAP (ABC)';
'GAP (ABC) <-> BPG (ABC)';
'BPG (ABC) <-> PGA (ABC)';
'PGA (ABC) <-> PEP (ABC)';
'PEP (ABC) <-> PYR (ABC)';
'm6PG (ABCDEF) <-> Ru5P (BCDEF) + CO2 (A)';
'Ru5P (ABCDE) <-> R5P (ABCDE)';
'Ru5P (ABCDE) <-> X5P (ABCDE)';
'X5P (ABCDE) + E4P (abcd) <-> GAP (CDE) + F6P (ABabcd)';
'X5P (ABCDE) + R5P (abcde) <-> S7P (ABabcde) + GAP (CDE)';
'GAP (ABC) + S7P (abcdefg) <-> E4P (defg) + F6P (abcABC)';
'DHAP (CBA) + E4P (DEFG) <-> SBP (ABCDEFG)';
'PYR (ABC) -> AcCoA (BC) + CO2 (A)';
'CitICit (ABCDEF) <-> OGA (ABCEF) + CO2 (D)';
'OGA (ABCDE) <-> SuccCoA (BCDE) + CO2 (A)';
'SuccCoA (ABCD) <-> Succ (ABCD)';
'Succ (ABCD) <-> Fum (ABCD)';
'Fum (ABCD) <-> MAL (ABCD)';
'MAL (ABCD) <-> OAA (ABCD)';
'CitICit (ABCDEF) -> Glx (AB) + Succ (DCEF)';
'Glx (AB) + AcCoA (ab) -> MAL (ABba)';
'CO2 (A) -> CO2.ex (A)';
'PGA (ABC) -> PSer (ABC)';
'PSer (ABC) -> Ser (ABC)';
'Ser (ABC) -> Ser.ex (ABC)';
'Ser (ABC) -> PYR (ABC)';
'Ser (ABC) <-> Gly (AB) + FTHF (C)';
'Gly (AB) <-> FTHF (B) + CO2 (A)';
'FTHF (A) -> FTHF.ex (A)';
'Gly (AB) -> Gly.ex (AB)';
'Ser (ABC) -> Cys (ABC)';
'Cys (ABC) -> Cys.ex (ABC)';
'PYR (ABC) -> Ala (ABC)';
'PYR (ABC) + PYR (abc) -> AKV (ABbcC) + CO2 (a)';
'AKV (ABbcC) -> Val (ABbcC)';
'AKV (ABCDE) + AcCoA (ab) -> IPPM (AabBCDE)';
'IPPM (AabBCDE) -> Leu (abBCDE) + CO2 (A)';
'Val (ABCDE) -> Val.ex (ABCDE)';
'Ala (ABC) -> Ala.ex (ABC)';
'Leu (ABCDEF) -> Leu.ex (ABCDEF)';
'PEP (ABC) + E4P (abcd) -> SKM (ABCabcd)';
'PEP (abc) + SKM (ABCDEFG) -> Chor (abcABCDEFG)';
'Chor (ABCDEFGHIJ) -> PHPYR (ABCEJIHGF) + CO2 (D)';
'PHPYR (ABCDEFGHI) -> Phe (ABCDEFGHI)';
'Chor (ABCDEFGHIJ) -> Tyr (ABCEFGHIJ) + CO2 (D)';
'Phe (ABCEFGHIJ) -> Phe.ex (ABCEFGHIJ)';
'Tyr (ABCEFGHIJ) -> Tyr.ex (ABCEFGHIJ)';
'Chor (abcABCDEFG) + PRPP (HIJKL) + Ser (MNO) -> Trp (MNOIHBCDEFG) + GAP (JKL) + PYR (abc) + CO2 (A)';
'Trp (ABCDEFGHIJK) -> Trp.ex (ABCDEFGHIJK)';
'OGA (ABCDE) -> Glu (ABCDE)';
'Glu (ABCDE) -> Pro (ABCDE)';
'Glu (ABCDE) -> Gln (ABCDE)';
'Glu (ABCDE) + AcCoA (ab) -> AcGlu (abABCDE)';
'AcGlu (abABCDE) -> AcCoA (ab) + ORN (ABCDE)';
'ORN (ABCDE) + CO2 (f) -> CITRL (ABCDEf)';
'CITRL (ABCDEf) -> Arg (ABCDEf)';
'Glu (ABCDE) -> Glu.ex (ABCDE)';
'Pro (ABCDE) -> Pro.ex (ABCDE)';
'Gln (ABCDE) -> Gln.ex (ABCDE)';
'Arg (ABCDEF) -> Arg.ex (ABCDEF)';
'Thr (ABCD) -> Gly (AB) + AcCoA (CD)';
'OAA (ABCD) -> Asp (ABCD)';
'Asp (ABCD) -> HomoSer (ABCD)';
'HomoSer (ABCD) -> Thr (ABCD)';
'Asp (ABCD) -> Asn (ABCD)';
'Thr (ABCD) + PYR (abc) -> Ile (ABbCDc) + CO2 (a)';
'HomoSer (ABCD) + FTHF (a) -> Met (ABCDa)';
'Asp (ABCD) + PYR (abc) -> Lys (ABCDcb) + CO2 (a)';
'Lys (ABCDEF) -> Lys.ex (ABCDEF)';
'Asp (ABCD) -> Asp.ex (ABCD)';
'Thr (ABCD) -> Thr.ex (ABCD)';
'Ile (ABCDEF) -> Ile.ex (ABCDEF)';
'Asn (ABCD) -> Asn.ex (ABCD)';
'Met (ABCDE) -> Met.ex (ABCDE)';
'Asp (ABCD) + CO2 (E) -> CBASP (ABCDE)';
'CBASP (ABCDE) -> OROT (ABCDE)';
'OROT (ABCDE) + PRPP (abcde) -> UMP (abcdeBCDE) + CO2 (A)';
'UMP (abcdeBCDE) + FTHF (f) -> dTTP (abcdeBCDEf)';
'UMP (abcdeBCDE) -> UMP.ex (abcdeBCDE)';
'dTTP (abcdeBCDEf) -> dTTP.ex (abcdeBCDEf)';
'G6P (ABCDEF) -> G6P.ex (ABCDEF)';
'F6P (ABCDEF) -> GAM6P (ABCDEF)';
'GAM6P (ABCDEF) + AcCoA (ab) -> ACGAM1P (ABCDEFab)';
'ACGAM1P (ABCDEFab) -> ACGAM1P.ex (ABCDEFab)';
'DHAP (ABC) -> GLYC3P (ABC)';
'GLYC3P (ABC) -> GLYC3P.ex (ABC)';
'R5P (ABCDE) -> PRPP (ABCDE)';
'PRPP (ABCDE) + ATP (abcdefghij) -> His (EDCBAj) + AICAR (abcdefghi)';
'His (ABCDEF) -> His.ex (ABCDEF)';
'PRPP (ABCDE) + Gly (gh) + FTHF (f) + CO2 (i) -> AICAR (ABCDEfghi)';
'AICAR (abcdefghi) + FTHF (J) -> ATP (abcdefghiJ)';
'ATP (ABCDEFGHIJ) -> ATP.ex (ABCDEFGHIJ)';
'AcCoA (AB) -> AcCoA.ex (AB)';
'Succ (ABCD) -> Succ.ex (ABCD)';
'PEP (ABC) + CO2 (D) -> OAA (ABCD)';
'MAL (ABCD) -> PYR (ABC) + CO2 (D)';
'OAA (ABCD) -> PYR (ABC) + CO2 (D)';
'G6P (ABCDEF) -> m6PG (ABCDEF)';
'm6PG (ABCDEF) -> PYR (ABC) + GAP (DEF)';
'SBP (ABCDEFG) -> S7P (ABCDEFG)';
'OAA (ABCD) + AcCoA (ab) -> CitICit (DCBAba)';
};

%% ---------- MS data ----------
msSpec = {
'CO20: CO20 @ 1';
'CO2: CO2 @ 1';
'GLC12: GLC12 @ 1 2 3 4 5 6';
'GLC: GLC @ 1 2 3 4 5 6';
'G6P: G6P @ 1 2 3 4 5 6';
'F6P: F6P @ 1 2 3 4 5 6';
'FBP: FBP @ 1 2 3 4 5 6';
'DHAP: DHAP @ 1 2 3';
'GAP: GAP @ 1 2 3';
'BPG: BPG @ 1 2 3';
'PGA: PGA @ 1 2 3';
'PEP: PEP @ 1 2 3';
'PYR: PYR @ 1 2 3';
'm6PG: m6PG @ 1 2 3 4 5 6';
'Ru5P: Ru5P @ 1 2 3 4 5';
'R5P: R5P @ 1 2 3 4 5';
'X5P: X5P @ 1 2 3 4 5';
'E4P: E4P @ 1 2 3 4';
'S7P: S7P @ 1 2 3 4 5 6 7';
'SBP: SBP @ 1 2 3 4 5 6 7';
'AcCoA: AcCoA @ 1 2';
'CitICit: CitICit @ 1 2 3 4 5 6';
'OGA: OGA @ 1 2 3 4 5';
'SuccCoA: SuccCoA @ 1 2 3 4';
'Succ: Succ @ 1 2 3 4';
'Fum: Fum @ 1 2 3 4';
'MAL: MAL @ 1 2 3 4';
'OAA: OAA @ 1 2 3 4';
'Glx: Glx @ 1 2';
'PSer: PSer @ 1 2 3';
'Ser: Ser @ 1 2 3';
'Gly: Gly @ 1 2';
'FTHF: FTHF @ 1';
'Cys: Cys @ 1 2 3';
'Ala: Ala @ 1 2 3';
'AKV: AKV @ 1 2 3 4 5';
'Val: Val @ 1 2 3 4 5';
'IPPM: IPPM @ 1 2 3 4 5 6 7';
'Leu: Leu @ 1 2 3 4 5 6';
'SKM: SKM @ 1 2 3 4 5 6 7';
'Chor: Chor @ 1 2 3 4 5 6 7 8 9 10';
'PHPYR: PHPYR @ 1 2 3 4 5 6 7 8 9';
'Phe: Phe @ 1 2 3 4 5 6 7 8 9';
'Tyr: Tyr @ 1 2 3 4 5 6 7 8 9';
'PRPP: PRPP @ 1 2 3 4 5';
'Trp: Trp @ 1 2 3 4 5 6 7 8 9 10 11';
'Glu: Glu @ 1 2 3 4 5';
'Pro: Pro @ 1 2 3 4 5';
'Gln: Gln @ 1 2 3 4 5';
'AcGlu: AcGlu @ 1 2 3 4 5 6 7';
'ORN: ORN @ 1 2 3 4 5';
'CITRL: CITRL @ 1 2 3 4 5 6';
'Arg: Arg @ 1 2 3 4 5 6';
'Thr: Thr @ 1 2 3 4';
'Asp: Asp @ 1 2 3 4';
'HomoSer: HomoSer @ 1 2 3 4';
'Asn: Asn @ 1 2 3 4';
'Ile: Ile @ 1 2 3 4 5 6';
'Met: Met @ 1 2 3 4 5';
'Lys: Lys @ 1 2 3 4 5 6';
'CBASP: CBASP @ 1 2 3 4 5';
'OROT: OROT @ 1 2 3 4 5';
'UMP: UMP @ 1 2 3 4 5 6 7 8 9';
'dTTP: dTTP @ 1 2 3 4 5 6 7 8 9 10';
'GAM6P: GAM6P @ 1 2 3 4 5 6';
'ACGAM1P: ACGAM1P @ 1 2 3 4 5 6 7 8';
'GLYC3P: GLYC3P @ 1 2 3';
'ATP: ATP @ 1 2 3 4 5 6 7 8 9 10';
'His: His @ 1 2 3 4 5 6';
'AICAR: AICAR @ 1 2 3 4 5 6 7 8 9';
};

%% ---------- Flux + tracers ----------
S = load('fluxfile.mat','fluxset'); % safe struct load inside a function/nested workspace
if ~isfield(S,'fluxset')
    error('fluxfile.mat does not contain variable ''fluxset''.');
end
Flux = S.fluxset;                 % [numSets x nFlux]

Tracers = { ...
    '[U-13C6]-glucose: GLC12 @ 1 2 3 4 5 6'; ...
    '[U-13C6]-glucose: GLC12 @ 1 2 3 4 5 6'; ...
    '[1-13C]-glucose : GLC12 @ 1'; ...
    '[2-13C]-glucose : GLC12 @ 2'; ...
    '[1,2-13C]-glucose : GLC12 @ 1 2'; ...
    '[6-13C]-glucose : GLC12 @ 6'; ...
    '[1,6-13C]-glucose : GLC12 @ 1 6'; ...
    '[3-13C]-glucose : GLC12 @ 3'; ...
    '[5-13C]-glucose : GLC12 @ 5'; ...
    '[4-13C]-glucose : GLC12 @ 4'; ...
    '[5,6-13C]-glucose : GLC12 @ 5 6'; ...
    '[1,2,3-13C]-glucose : GLC12 @ 1 2 3'; ...
    '[3,4-13C]-glucose : GLC12 @ 3 4' ...
};
frac = [1 0.5 1 1 1 1 1 1 1 1 1 1 1];

%% ---------- Ensure PROCESS-BASED pool ----------
p = gcp('nocreate');
if isempty(p)
    %parpool('local');
    delete(gcp('nocreate'));

    c = parcluster('Processes');   % matches the error profile
    c.NumWorkers = 2;              % start small
    
    % Put job storage somewhere writable and local
    c.JobStorageLocation = fullfile(tempdir, 'matlab_parpool_jobs');
    
    if ~exist(c.JobStorageLocation, 'dir')
        mkdir(c.JobStorageLocation);
    end
    
    parpool(c, 2);
elseif isa(p,'parallel.ThreadPool')
    delete(p); parpool('local');
end

%% ---------- Build worker-local model once per worker ----------
baseSpec = struct('Tracers',{Tracers}, 'frac', frac, 'msSpec',{msSpec}, ...
                  'reactionSpec',{reactionSpec});

% Preflight on client so errors are clear before going parallel
try
    tmpModel= buildWorkerModel(baseSpec); %#ok<NASGU>
catch ME
    error('Model build failed on client (before parallel): %s', ME.message);
end

modC = parallel.pool.Constant(@() buildWorkerModel(baseSpec), @cleanupWorkerModel);

%% ---------- Chunked parallel simulate over flux sets ----------
%numSets = size(Flux,1);
numSets = size(Flux,1);
if ~isempty(numSetsOverride)
    numSets = min(numSetsOverride, numSets);
end

chunkSize = 10000;                      % save/append every 1000 sims
outTxt    = 'raw_vals.txt';     % append after first chunk
chunkMat  = 'sim_raw_chunk.mat';       % overwritten every chunk

tStartAll = tic;
tParAll = 0; tSaveAll = 0; tWriteAll = 0;

fprintf('Running %d flux sets in chunks of %d ...\n', numSets, chunkSize);

for bStart = 1:chunkSize:numSets
    bEnd = min(bStart + chunkSize - 1, numSets);
    idx  = bStart:bEnd;
    nB   = numel(idx);

    % Per-chunk containers
    rawS_chunk    = cell(nB,1);
    timings_chunk = zeros(nB,1);

    % --- timing parfor for this chunk ---
    tPar = tic;
    parfor ii = 1:nB
        j = idx(ii)
        mloc = modC.Value              % worker-local model copy
        mloc.rates.flx.val = Flux(j,:); % set flux vector
        t0 = tic;
        Sj = simulate(mloc);
        timings_chunk(ii) = toc(t0);
        rawS_chunk{ii}    = Sj;
    end
    parSeconds = toc(tPar);
    tParAll = tParAll + parSeconds;
    fprintf('Chunk %d–%d: parfor %.3f s\n', bStart, bEnd, parSeconds);

    % --- save/overwrite .mat for this chunk (no I/O in parfor) ---
    tSave = tic;
    Flux_idx  = Flux(idx,:);   %#ok<NASGU>
    Tracers_b = Tracers;       %#ok<NASGU>
    msSpec_b  = msSpec;        %#ok<NASGU>
    timings   = timings_chunk; %#ok<NASGU>
    rawS      = rawS_chunk;    %#ok<NASGU>
    save(chunkMat,'rawS','Flux_idx','Tracers_b','msSpec_b','timings','idx','-v7.3');
    saveSeconds = toc(tSave);
    tSaveAll = tSaveAll + saveSeconds;
    fprintf('Chunk %d–%d: saved %s in %.3f s (overwritten)\n', bStart, bEnd, chunkMat, saveSeconds);

    % --- write/append named .val entries for this chunk ---
    mode = 'a'; if bStart == 1, mode = 'w'; end
    fid = fopen(outTxt, mode);
    if fid == -1, error('Could not open %s', outTxt); end
    c = onCleanup(@() fclose(fid));

    tWrite = tic;
    for ii = 1:nB
        j = idx(ii);
        fprintf(fid, '==== FluxSet %d ====\n', j);
        Sblk = rawS_chunk{ii};
        try
            if ~iscell(Sblk)
                SA = try_struct_array(Sblk);
                if ~isempty(SA)
                    for k = 1:numel(SA)
                        id = safe_get(SA(k),'id');
                        v  = safe_get(SA(k),'val');
                        label = label_from_id(id, k);
                        write_named_vals(fid, label, v);
                    end
                    fprintf(fid,'\n');
                    continue
                end
            end
            if iscell(Sblk)
                for w = 1:numel(Sblk)
                    fprintf(fid,'-- exp %d --\n', w);
                    s = Sblk{w};
                    SA = try_struct_array(s);
                    if ~isempty(SA)
                        for k = 1:numel(SA)
                            id = safe_get(SA(k),'id');
                            v  = safe_get(SA(k),'val');
                            label = label_from_id(id, k);
                            write_named_vals(fid, label, v);
                        end
                        continue
                    end
                    hasId  = (isstruct(s) && isfield(s,'id'))  || (isobject(s) && isprop(s,'id'));
                    hasVal = (isstruct(s) && isfield(s,'val')) || (isobject(s) && isprop(s,'val'));
                    if hasVal
                        vals = s.val; ids = [];
                        if hasId, ids = s.id; end
                        if iscell(vals)
                            for iii = 1:numel(vals)
                                lbl = label_from_id(index_cell(ids,iii), iii);
                                write_named_vals(fid, lbl, vals{iii});
                            end
                        else
                            lbl = label_from_id(ids, 1);
                            write_named_vals(fid, lbl, vals);
                        end
                        continue
                    end
                    dump = strtrim(evalc('disp(s)'));
                    fprintf(fid,'%s\n', dump);
                end
                fprintf(fid,'\n');
                continue
            end
            dump = strtrim(evalc('disp(Sblk)'));
            fprintf(fid,'%s\n\n', dump);
        catch ME
            fprintf(fid, '<<Error extracting vals for FluxSet %d: %s>>\n\n', j, ME.message);
        end
    end
    writeSeconds = toc(tWrite);
    tWriteAll = tWriteAll + writeSeconds;
    fprintf('Chunk %d–%d: appended to %s in %.3f s\n', bStart, bEnd, outTxt, writeSeconds);
end

%% ---------- Done ----------
fprintf('Parallel total: %.3f s | Save total: %.3f s | Write total: %.3f s\n', tParAll, tSaveAll, tWriteAll);
fprintf('Total elapsed time: %.2f s\n', toc(tStartAll));
end % main function

% ======================================================================
%                            SUBFUNCTIONS
% ======================================================================
function m = buildWorkerModel(spec)
% Per-worker model/experiment construction (no I/O)
    r = reaction(spec.reactionSpec);
    m = model(r);
    m.options.sim_ss = true;
    m.options.sim_na = false;

    % Symmetries
    m.mets{'Succ'}.sym = list('rotate180', atommap('1:4 2:3 3:2 4:1'));
    m.mets{'Fum'}.sym = list('rotate180', atommap('1:4 2:3 3:2 4:1'));

    d = msdata(spec.msSpec);

    % Build experiments (tracers)
    numT = numel(spec.Tracers);
    expts = repmat(experiment(), 1, numT);
    for w = 1:numT
        t = tracer(spec.Tracers{w});  % <<-- IMPORTANT: {} not ()
        t.frac = spec.frac(w);
        x = experiment(t);
        x.data_ms = d;
        expts(w) = x;
    end
    m.expts = expts;
end

function cleanupWorkerModel(~)
% placeholder (nothing to do); kept for symmetry & future resources
end

function SA = try_struct_array(x)
% Convert value/object/struct array to plain struct array, or [] on failure.
    SA = [];
    try
        SA = struct(x);
        return
    catch
    end
    try
        if isobject(x)
            mc = metaclass(x);
            hasVal = any(strcmp({mc.PropertyList.Name}, 'val'));
            hasId  = any(strcmp({mc.PropertyList.Name}, 'id'));
            if hasVal || hasId
                n = numel(x); SA = repmat(struct('id',[],'val',[]), n, 1);
                for i=1:n
                    SA(i).id  = try_get(x(i),'id');
                    SA(i).val = try_get(x(i),'val');
                end
            end
        end
    catch
        SA = [];
    end
end

function v = try_get(s, fld)
    v = [];
    try
        if isstruct(s) && isfield(s,fld), v = s.(fld); return; end
        if isobject(s)
            mc = metaclass(s);
            if any(strcmp({mc.PropertyList.Name}, fld))
                v = s.(fld); return
            end
        end
    catch
    end
end

function v = safe_get(s, fld)
    v = try_get(s, fld);
end

function lbl = label_from_id(id, k)
    if nargin < 2, k = 1; end
    if iscell(id), id = index_cell(id,1); end
    if isstring(id), id = char(id); end
    if ischar(id) && ~isempty(id)
        lbl = string(strtrim(id));
        return
    end
    try
        txt = strtrim(evalc('disp(id)'));
        if ~isempty(txt)
            lbl = string(txt); return
        end
    catch
    end
    lbl = "elem " + string(k);
end

function x = index_cell(c, i)
    if iscell(c) && i >= 1 && i <= numel(c)
        x = c{i};
    else
        x = [];
    end
end

function write_named_vals(fid, label, v)
    if isstring(label), label = char(label); end
    if ~ischar(label),  label = '<id>';      end
    fprintf(fid, '%s:', label);

    if iscell(v)
        nums = [];
        for ii = 1:numel(v)
            if isnumeric(v{ii}), nums = [nums, v{ii}(:).']; end %#ok<AGROW>
        end
        if ~isempty(nums)
            fprintf(fid, ' %s\n', sprintf('%.6f ', nums));
            return
        end
        raw = strtrim(evalc('disp(v)'));
        fprintf(fid, ' %s\n', raw);
        return
    end

    if isnumeric(v)
        fprintf(fid, ' %s\n', sprintf('%.6f ', v(:).'));
        return
    end

    if isstring(v) || ischar(v)
        fprintf(fid, ' %s\n', string(v));
        return
    end

    raw = strtrim(evalc('disp(v)'));
    fprintf(fid, ' %s\n', raw);
end
