# Скрипт для инициализации conda в текущей сессии PowerShell
# Запустите этот скрипт в начале каждой сессии или добавьте в профиль PowerShell

# Добавляем conda в PATH текущей сессии
$condaScripts = "$env:LOCALAPPDATA\miniconda3\Scripts"
$condaLibrary = "$env:LOCALAPPDATA\miniconda3\Library\bin"

if (Test-Path $condaScripts) {
    if ($env:PATH -notlike "*$condaScripts*") {
        $env:PATH = "$condaScripts;$condaLibrary;$env:PATH"
        Write-Host "Conda добавлена в PATH текущей сессии" -ForegroundColor Green
    }
    
    # Инициализируем conda для PowerShell
    & "$condaScripts\conda.exe" init powershell --quiet
    
    Write-Host "Conda инициализирована. Перезапустите PowerShell или выполните:" -ForegroundColor Yellow
    Write-Host "  . `$PROFILE" -ForegroundColor Cyan
} else {
    Write-Host "Conda не найдена в $condaScripts" -ForegroundColor Red
}



