# Скрипт для добавления conda в PATH

$condaBasePath = Join-Path $env:LOCALAPPDATA "miniconda3"

if (Test-Path $condaBasePath) {
    Write-Host "Conda найдена в: $condaBasePath" -ForegroundColor Green
    
    $scriptsPath = Join-Path $condaBasePath "Scripts"
    $libraryPath = Join-Path $condaBasePath "Library\bin"
    
    $currentPath = [Environment]::GetEnvironmentVariable("Path", "User")
    $pathUpdated = $false
    
    if (-not $currentPath.Contains($scriptsPath)) {
        $currentPath = $currentPath + ";" + $scriptsPath
        $pathUpdated = $true
        Write-Host "Добавлен в PATH: $scriptsPath" -ForegroundColor Yellow
    }
    
    if (-not $currentPath.Contains($libraryPath)) {
        $currentPath = $currentPath + ";" + $libraryPath
        $pathUpdated = $true
        Write-Host "Добавлен в PATH: $libraryPath" -ForegroundColor Yellow
    }
    
    if ($pathUpdated) {
        [Environment]::SetEnvironmentVariable("Path", $currentPath, "User")
        Write-Host ""
        Write-Host "PATH успешно обновлен!" -ForegroundColor Green
        Write-Host "Перезапустите PowerShell для применения изменений" -ForegroundColor Cyan
    } else {
        Write-Host ""
        Write-Host "PATH уже содержит все необходимые пути" -ForegroundColor Green
    }
    
    $condaExe = Join-Path $scriptsPath "conda.exe"
    if (Test-Path $condaExe) {
        Write-Host ""
        Write-Host "Проверка conda:" -ForegroundColor Cyan
        & $condaExe --version
    }
} else {
    Write-Host "Conda не найдена в: $condaBasePath" -ForegroundColor Red
}
