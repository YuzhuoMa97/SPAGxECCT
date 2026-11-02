@echo off
setlocal
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem Optional CCT; run tutorial 04 first. Match metadata, phenotype and environment.
rem Matched total G is a toy sum of local dosages. Real admixed analyses require suitable SNP PCs.
"%SPAGXE_EXE%" --vcf examples/local/global.vcf --pheno examples/local/samples.pheno --pheno-name Quantitative --covar-name X,E --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/global
if errorlevel 1 exit /b 1
rem Dynamic weights give one nominal weight to each valid ancestry and one to global.
"%SPAGXE_EXE%" --method SPAGxEcct --cct-global results/tutorials/global.Quantitative.SPAGxE --cct-local results/tutorials/local_text.Quantitative.SPAGxE_local_CCT --cct-phenotype Quantitative --envir-name E --cct-per-ancestry-equal --out results/tutorials/global_local
if errorlevel 1 exit /b 1
exit /b 0
