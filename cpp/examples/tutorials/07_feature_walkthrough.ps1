# Complete feature walkthrough. All analyses call the standalone executable.
# File-list preparation uses built-in PowerShell; no R or Python is needed.
param([string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$cppDirectory = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../..'))
$program = $env:SPAGXE_EXE
if (-not $program) { $program = $env:SPAGXE_BIN }
if (-not $program) {
    foreach ($relative in @('SPAGxE.exe', 'build/bin/SPAGxE.exe')) {
        $candidate = Join-Path $cppDirectory $relative
        if (Test-Path -LiteralPath $candidate -PathType Leaf) { $program = $candidate; break }
    }
}
if (-not $program -or -not (Test-Path -LiteralPath $program -PathType Leaf)) {
    throw 'Set SPAGXE_EXE to a built SPAGxE executable.'
}
$program = (Resolve-Path -LiteralPath $program).Path
if (-not $OutputDirectory) { $OutputDirectory = 'results/official_walkthrough' }
if (-not [IO.Path]::IsPathRooted($OutputDirectory)) {
    $OutputDirectory = Join-Path $cppDirectory $OutputDirectory
}
$outDirectory = [IO.Path]::GetFullPath($OutputDirectory)
if (Test-Path -LiteralPath $outDirectory) {
    if (-not (Test-Path -LiteralPath $outDirectory -PathType Container) -or
        @(Get-ChildItem -LiteralPath $outDirectory -Force).Count -gt 0) {
        throw "Use a new or empty output directory: $outDirectory"
    }
}
[IO.Directory]::CreateDirectory($outDirectory) | Out-Null
$utf8 = New-Object Text.UTF8Encoding($false)
function ConvertTo-NativeArgument([string]$Value) {
    # Windows native quoting preserves spaces, quotes and trailing backslashes.
    $escaped = [regex]::Replace($Value, '(\\*)"', '$1$1\\"')
    $escaped = [regex]::Replace($escaped, '(\\+)$', '$1$1')
    return '"' + $escaped + '"'
}
function Invoke-SPAGxE([string]$Example, [string[]]$Arguments) {
    $startInfo = New-Object Diagnostics.ProcessStartInfo
    $startInfo.FileName = $program
    $startInfo.Arguments = (($Arguments | ForEach-Object { ConvertTo-NativeArgument $_ }) -join ' ')
    $startInfo.WorkingDirectory = $cppDirectory
    $startInfo.UseShellExecute = $false
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.StandardOutputEncoding = $utf8
    $startInfo.StandardErrorEncoding = $utf8
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $startInfo
    try {
        if (-not $process.Start()) { throw "Could not start $Example" }
        $stdout = $process.StandardOutput.ReadToEndAsync()
        $stderr = $process.StandardError.ReadToEndAsync()
        $process.WaitForExit()
        [Console]::Out.Write($stdout.Result)
        [Console]::Error.Write($stderr.Result)
        if ($process.ExitCode -ne 0) { throw "$Example failed with exit code $($process.ExitCode)" }
    } finally { $process.Dispose() }
}


# UG01_help: Executable version, general help, all options, phenotype and output help
Write-Host 'UG01_help'
Invoke-SPAGxE -Example 'UG01_help' -Arguments @(
    '--version'
)
Invoke-SPAGxE -Example 'UG01_help' -Arguments @(
    '--help'
)
Invoke-SPAGxE -Example 'UG01_help' -Arguments @(
    '--help',
    'options'
)
Invoke-SPAGxE -Example 'UG01_help' -Arguments @(
    '--help',
    'phenotype'
)
Invoke-SPAGxE -Example 'UG01_help' -Arguments @(
    '--help',
    'output'
)

# UG02_base: SPAGxE native quantitative null, PLINK BED, residual saving
Write-Host 'UG02_base'
Invoke-SPAGxE -Example 'UG02_base' -Arguments @(
    '--method',
    'SPAGxE',
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--regression-model',
    'linear',
    '--save-resid',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'base')
)

# UG03_plus_table: SPAGxE+ shared-AF GRM correction using ID1 ID2 VALUE table
Write-Host 'UG03_plus_table'
Invoke-SPAGxE -Example 'UG03_plus_table' -Arguments @(
    '--method',
    'SPAGxE+',
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--grm-file',
    'examples/synthetic/synthetic.grm.tsv',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'plus_table')
)

# UG04_mix: SPAGxEmix individual-specific AF, ancestry PCs, embedded-ID BGEN ref-first
Write-Host 'UG04_mix'
Invoke-SPAGxE -Example 'UG04_mix' -Arguments @(
    '--method',
    'SPAGxEmix',
    '--bgen',
    'examples/synthetic/synthetic.bgen',
    'ref-first',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--pc-cols',
    'PC1,PC2',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'mix')
)

# UG05_local_text: SPAGxEmix_local paired G/H, four native families, across-ancestry CCT
Write-Host 'UG05_local_text'
Invoke-SPAGxE -Example 'UG05_local_text' -Arguments @(
    '--method',
    'SPAGxEmix_local',
    '--local-manifest',
    'examples/local/cohort.local.tsv',
    '--pheno',
    'examples/local/samples.pheno',
    '--pheno-name',
    'Quantitative,Binary,Time:Event,Ordinal',
    '--covar-name',
    'X,E',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'local_text')
)

# UG06_bgen_ref_last: Global BGEN last-listed REF orientation
Write-Host 'UG06_bgen_ref_last'
Invoke-SPAGxE -Example 'UG06_bgen_ref_last' -Arguments @(
    '--bgen',
    'examples/synthetic/synthetic.reflast.bgen',
    'ref-last',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'bgen_ref_last')
)

# UG33_bgen_external_ids: Global BGEN without embedded IIDs, companion sample resolution
Write-Host 'UG33_bgen_external_ids'
Invoke-SPAGxE -Example 'UG33_bgen_external_ids' -Arguments @(
    '--bgen',
    'examples/synthetic/synthetic.noids.bgen',
    'ref-first',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'bgen_external_ids')
)

# UG07_vcf: VCF global reader
Write-Host 'UG07_vcf'
Invoke-SPAGxE -Example 'UG07_vcf' -Arguments @(
    '--vcf',
    'examples/synthetic/synthetic.vcf',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'vcf')
)

# UG08_bcf: BCF2 global reader
Write-Host 'UG08_bcf'
Invoke-SPAGxE -Example 'UG08_bcf' -Arguments @(
    '--bcf',
    'examples/synthetic/synthetic.bcf',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'bcf')
)

# UG09_pgen: PGEN/PVAR/PSAM global reader with matched 1kg phenotype fixture
Write-Host 'UG09_pgen'
Invoke-SPAGxE -Example 'UG09_pgen' -Arguments @(
    '--pfile',
    'examples/1kg',
    '--pheno',
    'examples/1kg.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'MALE,PC1,PC2,PC3,PC4',
    '--envir-name',
    'MALE',
    '--out',
    (Join-Path $outDirectory 'pgen')
)

# UG10_plus_plink2: SPAGxE+ PLINK2 sparse GRM plus grm.id
Write-Host 'UG10_plus_plink2'
Invoke-SPAGxE -Example 'UG10_plus_plink2' -Arguments @(
    '--method',
    'SPAGxE+',
    '--pfile',
    'examples/1kg',
    '--pheno',
    'examples/1kg.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'MALE,PC1,PC2,PC3,PC4',
    '--envir-name',
    'MALE',
    '--sp-grm-plink2',
    'examples/1kg.grm.sp',
    '--out',
    (Join-Path $outDirectory 'plus_plink2')
)

# UG11_trait_models: Explicit native logistic, Breslow Cox and ordinal fitting
Write-Host 'UG11_trait_models'
Invoke-SPAGxE -Example 'UG11_trait_models' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Binary',
    '--regression-model',
    'logistic',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'binary')
)
Invoke-SPAGxE -Example 'UG11_trait_models' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Time:Event',
    '--regression-model',
    'cox',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'survival')
)
Invoke-SPAGxE -Example 'UG11_trait_models' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Ordinal',
    '--regression-model',
    'ordinal',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'ordinal')
)

# UG12_multiple: Multiple raw traits and two environments with automatic families
Write-Host 'UG12_multiple'
Invoke-SPAGxE -Example 'UG12_multiple' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative,Binary,Time:Event,Ordinal',
    '--regression-model',
    'auto',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E,X1',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'multiple')
)

# UG13_external_residuals: Several precomputed residual columns
Write-Host 'UG13_external_residuals'
Invoke-SPAGxE -Example 'UG13_external_residuals' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--resid-name',
    'R_quant,R_binary,R_survival,R_ordinal',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'external_residuals')
)

# UG14_reload_residual: Native saved-residual reload with separate environment table
Write-Host 'UG14_reload_residual'
Invoke-SPAGxE -Example 'UG14_reload_residual' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    (Join-Path $outDirectory 'base.null.resid'),
    '--resid-name',
    'Quantitative',
    '--covar',
    'examples/synthetic/synthetic.pheno',
    '--covar-name',
    'E',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'reloaded_residual')
)

# UG15_prepare_filters: Create IID and marker include/exclude lists
Write-Host 'UG15_prepare_filters'
$keepLines = [IO.File]::ReadAllLines((Join-Path $cppDirectory 'examples/synthetic/synthetic.pheno')) | Select-Object -First 401 | ForEach-Object { [regex]::Split($_.Trim(), '\s+')[0] }
[IO.File]::WriteAllText((Join-Path $outDirectory 'keep.txt'), (($keepLines -join "`n") + "`n"), $utf8)
[IO.File]::WriteAllText((Join-Path $outDirectory 'remove.txt'), 'IID
S0001
S0002
', $utf8)
[IO.File]::WriteAllText((Join-Path $outDirectory 'extract.txt'), 'rsOracle01
rsOracle02
rsOracle03
rsOracle04
', $utf8)
[IO.File]::WriteAllText((Join-Path $outDirectory 'exclude.txt'), 'rsOracle02
', $utf8)

# UG16_filter_qc: Combined keep/remove/extract/exclude/chromosome selection and all QC flags
Write-Host 'UG16_filter_qc'
Invoke-SPAGxE -Example 'UG16_filter_qc' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--keep',
    (Join-Path $outDirectory 'keep.txt'),
    '--remove',
    (Join-Path $outDirectory 'remove.txt'),
    '--extract',
    (Join-Path $outDirectory 'extract.txt'),
    '--exclude',
    (Join-Path $outDirectory 'exclude.txt'),
    '--chr',
    '1',
    '--geno',
    '0.1',
    '--maf',
    '0.01',
    '--mac',
    '10',
    '--hwe',
    '0.000001',
    '--hard-call-threshold',
    '0.1',
    '--out',
    (Join-Path $outDirectory 'filtered')
)

# UG17_local_existing_counts: Existing native packed G/H with known H and supplied residual
Write-Host 'UG17_local_existing_counts'
Invoke-SPAGxE -Example 'UG17_local_existing_counts' -Arguments @(
    '--method',
    'SPAGxEmix_local',
    '--local-prefix',
    'examples/local/binary_demo/local',
    '--pheno',
    'examples/local/binary_demo/samples.pheno',
    '--resid-name',
    'R',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'local_existing')
)

# UG18_local_bgen: Local known-H BGEN layout2 zlib/16-bit profile
Write-Host 'UG18_local_bgen'
Invoke-SPAGxE -Example 'UG18_local_bgen' -Arguments @(
    '--method',
    'SPAGxEmix_local',
    '--local-bgen-manifest',
    'examples/local/bgen_demo/bits16_zlib.localbgen.tsv',
    '--pheno',
    'examples/local/bgen_demo/samples.pheno',
    '--resid-name',
    'R',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'local_bgen')
)

# UG19_pack_local: Native C++ integer G/H packer, original BIM/FAM reuse, block/level/budget controls
Write-Host 'UG19_pack_local'
Invoke-SPAGxE -Example 'UG19_pack_local' -Arguments @(
    '--pack-local',
    '--local-manifest',
    'examples/local/binary_demo/cohort.local.tsv',
    '--local-bim',
    'examples/local/binary_demo/local.localbim',
    '--local-fam',
    'examples/local/binary_demo/local.localfam',
    '--local-block-variants',
    '2',
    '--local-zstd-level',
    '5',
    '--local-max-block-bytes',
    '268435456',
    '--out',
    (Join-Path $outDirectory 'local_counts')
)

# UG20_local_new_counts: Analyze newly generated native G/H triplet
Write-Host 'UG20_local_new_counts'
Invoke-SPAGxE -Example 'UG20_local_new_counts' -Arguments @(
    '--method',
    'SPAGxEmix_local',
    '--local-prefix',
    (Join-Path $outDirectory 'local_counts'),
    '--pheno',
    'examples/local/binary_demo/samples.pheno',
    '--resid-name',
    'R',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'local_new')
)

# UG21_repack_local: Repack existing native G/H triplet with new block size and default zstd level
Write-Host 'UG21_repack_local'
Invoke-SPAGxE -Example 'UG21_repack_local' -Arguments @(
    '--pack-local',
    '--local-prefix',
    (Join-Path $outDirectory 'local_counts'),
    '--local-block-variants',
    '1',
    '--local-zstd-level',
    '0',
    '--local-max-block-bytes',
    '268435456',
    '--out',
    (Join-Path $outDirectory 'local_repacked')
)

# UG22_cct_global_input: Global scan matching local fixture variant catalog for CCT
Write-Host 'UG22_cct_global_input'
Invoke-SPAGxE -Example 'UG22_cct_global_input' -Arguments @(
    '--vcf',
    'examples/local/global.vcf',
    '--pheno',
    'examples/local/samples.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X,E',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'cct_global')
)

# UG23_cct_default: Global/local CCT with default equal group weights
Write-Host 'UG23_cct_default'
Invoke-SPAGxE -Example 'UG23_cct_default' -Arguments @(
    '--method',
    'SPAGxEcct',
    '--cct-global',
    (Join-Path $outDirectory 'cct_global.Quantitative.SPAGxE'),
    '--cct-local',
    (Join-Path $outDirectory 'local_text.Quantitative.SPAGxE_local_CCT'),
    '--cct-phenotype',
    'Quantitative',
    '--envir-name',
    'E',
    '--out',
    (Join-Path $outDirectory 'cct_default')
)

# UG24_cct_weights: Global/local CCT custom positive normalized group weights
Write-Host 'UG24_cct_weights'
Invoke-SPAGxE -Example 'UG24_cct_weights' -Arguments @(
    '--method',
    'SPAGxEcct',
    '--cct-global',
    (Join-Path $outDirectory 'cct_global.Quantitative.SPAGxE'),
    '--cct-local',
    (Join-Path $outDirectory 'local_text.Quantitative.SPAGxE_local_CCT'),
    '--cct-phenotype',
    'Quantitative',
    '--envir-name',
    'E',
    '--cct-weights',
    '0.6,0.4',
    '--out',
    (Join-Path $outDirectory 'cct_weighted')
)

# UG25_cct_components: Global/local CCT equal nominal weight per valid component
Write-Host 'UG25_cct_components'
Invoke-SPAGxE -Example 'UG25_cct_components' -Arguments @(
    '--method',
    'SPAGxEcct',
    '--cct-global',
    (Join-Path $outDirectory 'cct_global.Quantitative.SPAGxE'),
    '--cct-local',
    (Join-Path $outDirectory 'local_text.Quantitative.SPAGxE_local_CCT'),
    '--cct-phenotype',
    'Quantitative',
    '--envir-name',
    'E',
    '--cct-per-ancestry-equal',
    '--out',
    (Join-Path $outDirectory 'cct_equal_components')
)

# UG26_af_export: Individual-AF export and gzip codec with default compression level
Write-Host 'UG26_af_export'
Invoke-SPAGxE -Example 'UG26_af_export' -Arguments @(
    '--cal-af-coef',
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pc-cols',
    'PC1,PC2',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--hwe',
    '0',
    '--threads',
    '2',
    '--compression',
    'gz',
    '--out',
    (Join-Path $outDirectory 'af_models')
)

# UG27_int: Native inverse-normal transform of one selected quantitative phenotype
Write-Host 'UG27_int'
Invoke-SPAGxE -Example 'UG27_int' -Arguments @(
    '--int-pheno',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--out',
    (Join-Path $outDirectory 'int_quantitative')
)

# UG28_separate_covariates: Analyze transformed phenotype with separately supplied original covariates
Write-Host 'UG28_separate_covariates'
Invoke-SPAGxE -Example 'UG28_separate_covariates' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    (Join-Path $outDirectory 'int_quantitative.txt'),
    '--pheno-name',
    'Quantitative',
    '--covar',
    'examples/synthetic/synthetic.pheno',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'int_analysis')
)

# UG29_gzip: Gzip association output plus saved residual output
Write-Host 'UG29_gzip'
Invoke-SPAGxE -Example 'UG29_gzip' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--save-resid',
    '--compression',
    'gz',
    '--compression-level',
    '6',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'gzip')
)

# UG30_zstd: Zstd association output and explicit level
Write-Host 'UG30_zstd'
Invoke-SPAGxE -Example 'UG30_zstd' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--compression',
    'zst',
    '--compression-level',
    '3',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'zstd')
)

# UG31_compressed_residual: Read gzip saved residual directly and separate environment table
Write-Host 'UG31_compressed_residual'
Invoke-SPAGxE -Example 'UG31_compressed_residual' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    (Join-Path $outDirectory 'gzip.null.resid.gz'),
    '--resid-name',
    'Quantitative',
    '--covar',
    'examples/synthetic/synthetic.pheno',
    '--covar-name',
    'E',
    '--envir-name',
    'E',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'compressed_residual')
)

# UG32_resources: Threads/chunks/seed compatibility plus explicit branch and SPA thresholds
Write-Host 'UG32_resources'
Invoke-SPAGxE -Example 'UG32_resources' -Arguments @(
    '--bfile',
    'examples/synthetic/synthetic',
    '--pheno',
    'examples/synthetic/synthetic.pheno',
    '--pheno-name',
    'Quantitative',
    '--covar-name',
    'X1,E,PC1,PC2',
    '--envir-name',
    'E',
    '--threads',
    '4',
    '--chunk-ksnp',
    '0.5',
    '--seed',
    '20261029',
    '--spagxe-marginal-cutoff',
    '0.001',
    '--spa-z-threshold',
    '2',
    '--mac',
    '0',
    '--maf',
    '0',
    '--geno',
    '1',
    '--out',
    (Join-Path $outDirectory 'resources')
)
Write-Host "Feature walkthrough completed: $outDirectory"
