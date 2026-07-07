$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoDir = Resolve-Path (Join-Path $ScriptDir "..\..")
$OutDir = Join-Path $RepoDir "03_result\Plasma_Protein\Sankey_color_optimized"
$Report = Join-Path $OutDir "report.txt"
$Rscript = "D:\R\R-4.5.1\bin\x64\Rscript.exe"
$Pipeline = Join-Path $ScriptDir "abundance_sankey_color_optimized.R"

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

"===== Pipeline run started: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') =====" | Tee-Object -FilePath $Report -Append
"Rscript: $Rscript" | Tee-Object -FilePath $Report -Append
"Pipeline: $Pipeline" | Tee-Object -FilePath $Report -Append

if (-not (Test-Path -LiteralPath $Rscript)) {
    "ERROR: Rscript not found at $Rscript" | Tee-Object -FilePath $Report -Append
    exit 1
}

if (-not (Test-Path -LiteralPath $Pipeline)) {
    "ERROR: Pipeline script not found at $Pipeline" | Tee-Object -FilePath $Report -Append
    exit 1
}

$command = "`"$Rscript`" `"$Pipeline`""
cmd.exe /c "$command 2>&1" | Tee-Object -FilePath $Report -Append
$ExitCode = $LASTEXITCODE

"===== Pipeline run finished: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss'); exit code: $ExitCode =====" | Tee-Object -FilePath $Report -Append
exit $ExitCode
