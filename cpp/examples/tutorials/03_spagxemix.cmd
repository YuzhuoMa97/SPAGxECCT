@echo off
setlocal
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem SPAGxEmix: individual AF from genotype PCs.
rem The same PCs enter AF and null models. BGEN allele orientation is explicit.
"%SPAGXE_EXE%" --method SPAGxEmix --bgen examples/synthetic/synthetic.bgen ref-first --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X1,E,PC1,PC2 --envir-name E --pc-cols PC1,PC2 --threads 4 --out results/tutorials/mix
if errorlevel 1 exit /b 1
exit /b 0
