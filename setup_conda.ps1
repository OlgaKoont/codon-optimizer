# Скрипт для настройки conda окружения kodon_opt
# Запустите этот скрипт после установки Miniconda/Anaconda

Write-Host "Проверка установки conda..." -ForegroundColor Yellow

# Проверяем различные возможные пути установки conda
$condaPaths = @(
    "$env:USERPROFILE\Miniconda3\Scripts\conda.exe",
    "$env:USERPROFILE\Anaconda3\Scripts\conda.exe",
    "$env:LOCALAPPDATA\miniconda3\Scripts\conda.exe",
    "$env:LOCALAPPDATA\anaconda3\Scripts\conda.exe",
    "$env:ProgramData\Miniconda3\Scripts\conda.exe",
    "$env:ProgramData\Anaconda3\Scripts\conda.exe"
)

$condaExe = $null
foreach ($path in $condaPaths) {
    if (Test-Path $path) {
        $condaExe = $path
        Write-Host "Найден conda: $condaExe" -ForegroundColor Green
        break
    }
}

if (-not $condaExe) {
    Write-Host "Conda не найдена. Пожалуйста:" -ForegroundColor Red
    Write-Host "1. Скачайте Miniconda с https://docs.conda.io/en/latest/miniconda.html" -ForegroundColor Yellow
    Write-Host "2. Установите Miniconda" -ForegroundColor Yellow
    Write-Host "3. Перезапустите PowerShell" -ForegroundColor Yellow
    Write-Host "4. Запустите этот скрипт снова" -ForegroundColor Yellow
    exit 1
}

# Инициализация conda для PowerShell
Write-Host "Инициализация conda для PowerShell..." -ForegroundColor Yellow
& $condaExe init powershell

Write-Host "`nСоздание окружения kodon_opt..." -ForegroundColor Yellow
& $condaExe create -n kodon_opt python=3.8 -y

Write-Host "`nАктивация окружения..." -ForegroundColor Yellow
Write-Host "Выполните следующую команду для активации:" -ForegroundColor Cyan
Write-Host "conda activate kodon_opt" -ForegroundColor Green

Write-Host "`nПосле активации выполните:" -ForegroundColor Cyan
Write-Host "pip install -r requirements.txt" -ForegroundColor Green
Write-Host "pip install -e ." -ForegroundColor Green

Write-Host "`nГотово! Перезапустите PowerShell для применения изменений." -ForegroundColor Green



