# Script to setup PowerShell profile for conda auto-initialization

$profilePath = $PROFILE

Write-Host "PowerShell profile path: $profilePath" -ForegroundColor Cyan

# Create profile if it doesn't exist
if (!(Test-Path $profilePath)) {
    $profileDir = Split-Path $profilePath -Parent
    if (!(Test-Path $profileDir)) {
        New-Item -Path $profileDir -ItemType Directory -Force | Out-Null
    }
    New-Item -Path $profilePath -ItemType File -Force | Out-Null
    Write-Host "Profile created" -ForegroundColor Green
}

# Check if conda initialization already exists
$profileContent = Get-Content $profilePath -ErrorAction SilentlyContinue
if ($profileContent -like "*miniconda3*") {
    Write-Host "Conda already configured in profile" -ForegroundColor Yellow
} else {
    # Add conda initialization
    $condaInit = @"

# Initialize conda
`$condaPath = "`$env:LOCALAPPDATA\miniconda3"
if (Test-Path "`$condaPath\Scripts\conda.exe") {
    `$env:PATH = "`$condaPath\Scripts;`$condaPath\Library\bin;`$env:PATH"
    (& "`$condaPath\Scripts\conda.exe" "shell.powershell" "hook") | Out-String | Invoke-Expression
}
"@
    
    Add-Content $profilePath $condaInit
    Write-Host "Conda initialization added to profile!" -ForegroundColor Green
    Write-Host "Restart PowerShell to apply changes" -ForegroundColor Yellow
}
