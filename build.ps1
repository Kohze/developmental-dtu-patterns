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
    & $pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
    if ($LASTEXITCODE -ne 0) { throw 'First LaTeX pass failed.' }
    & $bibtex manuscript
    if ($LASTEXITCODE -ne 0) { throw 'BibTeX failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
    if ($LASTEXITCODE -ne 0) { throw 'Second LaTeX pass failed.' }
    & $pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
    if ($LASTEXITCODE -ne 0) { throw 'Final LaTeX pass failed.' }
}
finally {
    Pop-Location
}

$patterns = 'undefined citations|undefined references|Overfull \\hbox|Fatal error'
$bad = Select-String -LiteralPath (Join-Path $paperDir 'manuscript.log') -Pattern $patterns -CaseSensitive:$false
if ($bad) { throw 'The manuscript built with unresolved log problems.' }

Get-Item -LiteralPath (Join-Path $paperDir 'manuscript.pdf')
