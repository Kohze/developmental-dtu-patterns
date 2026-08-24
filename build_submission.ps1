param([string]$TexBin)

$ErrorActionPreference = 'Stop'
$paperDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$localBin = Join-Path $paperDir 'tools\miktex\texmfs\install\miktex\bin\x64'
$selected = if ($TexBin) { (Resolve-Path -LiteralPath $TexBin).Path } else { $localBin }
$pdflatex = Join-Path $selected 'pdflatex.exe'
$bibtex = Join-Path $selected 'bibtex.exe'

if (-not (Test-Path -LiteralPath $pdflatex)) {
    $pdflatex = (Get-Command pdflatex.exe -ErrorAction Stop).Source
    $bibtex = (Get-Command bibtex.exe -ErrorAction Stop).Source
}

$env:PATH = "$(Split-Path -Parent $pdflatex);$env:PATH"
$env:SOURCE_DATE_EPOCH = '946684800'
$env:FORCE_SOURCE_DATE = '1'

Push-Location $paperDir
try {
    & $pdflatex -interaction=nonstopmode -halt-on-error submission_main.tex
    if ($LASTEXITCODE -ne 0) { throw 'First main-article LaTeX pass failed.' }
    & $bibtex submission_main
    if ($LASTEXITCODE -ne 0) { throw 'Main-article BibTeX pass failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error submission_main.tex
    if ($LASTEXITCODE -ne 0) { throw 'Second main-article LaTeX pass failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error submission_main.tex
    if ($LASTEXITCODE -ne 0) { throw 'Final main-article LaTeX pass failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error supplementary_figures.tex
    if ($LASTEXITCODE -ne 0) { throw 'Supplementary-figure LaTeX pass failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error supplementary_figures.tex
    if ($LASTEXITCODE -ne 0) { throw 'Final supplementary-figure LaTeX pass failed.' }
}
finally {
    Pop-Location
}

$patterns = 'undefined citations|undefined references|Citation.*undefined|Reference.*undefined|multiply defined|Overfull \\hbox|Fatal error|LaTeX Error'
$bad = Select-String -LiteralPath (Join-Path $paperDir 'submission_main.log'), (Join-Path $paperDir 'supplementary_figures.log') -Pattern $patterns -CaseSensitive:$false
if ($bad) { throw 'The submission files built with unresolved log problems.' }

Get-Item -LiteralPath (Join-Path $paperDir 'submission_main.pdf'), (Join-Path $paperDir 'supplementary_figures.pdf')
