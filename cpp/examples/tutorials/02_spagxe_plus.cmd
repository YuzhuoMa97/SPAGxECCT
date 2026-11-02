@echo off
setlocal
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem SPAGxE+: sparse-GRM covariance correction; SPA-only final P.
rem Use bundled public PGEN/phenotype/GRM inputs. Keep the companion .grm.id.
"%SPAGXE_EXE%" --method SPAGxE+ --pfile examples/1kg --pheno examples/1kg.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name MALE,PC1,PC2,PC3,PC4 --envir-name MALE --sp-grm-plink2 examples/1kg.grm.sp --threads 4 --out results/tutorials/plus
if errorlevel 1 exit /b 1
exit /b 0
