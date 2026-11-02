@echo off
setlocal
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem SPAGxEmix_local: all local examples are independently generated SPAGxE simulations.
rem Paired G/H tables; native four-family null models and Wald/CCT. Relaxed QC exposes tiny-example edge cases.
"%SPAGXE_EXE%" --method SPAGxEmix_local --local-manifest examples/local/cohort.local.tsv --pheno examples/local/samples.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --threads 4 --out results/tutorials/local_text
if errorlevel 1 exit /b 1
rem Packed hardcalls, H=0/1/2. External residuals have no prospective Wald leg.
"%SPAGXE_EXE%" --method SPAGxEmix_local --local-prefix examples/local/binary_demo/local --pheno examples/local/binary_demo/samples.pheno --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/local_binary
if errorlevel 1 exit /b 1
rem Known-H BGEN layout 2 per ancestry; use expected G from stored posteriors. No R runtime.
"%SPAGXE_EXE%" --method SPAGxEmix_local --local-bgen-manifest examples/local/bgen_demo/bits16_zlib.localbgen.tsv --pheno examples/local/samples.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --threads 4 --out results/tutorials/local_bgen
if errorlevel 1 exit /b 1
exit /b 0
