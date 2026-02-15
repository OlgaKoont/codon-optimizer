# Скрипт для активации окружения kodon_opt в текущей сессии
# Запустите: . .\activate_kodon_opt.ps1

# Обновляем PATH
$env:PATH = "$env:LOCALAPPDATA\miniconda3\Scripts;$env:LOCALAPPDATA\miniconda3\Library\bin;$env:PATH"

# Активируем окружение
& "$env:LOCALAPPDATA\miniconda3\Scripts\conda.exe" activate kodon_opt

# Создаем алиас для codon-optimize
function codon-optimize {
    & "$env:LOCALAPPDATA\miniconda3\envs\kodon_opt\python.exe" -m codon_optimizer.cli.main $args
}

Write-Host "Окружение kodon_opt активировано!" -ForegroundColor Green
Write-Host "Теперь можно использовать:" -ForegroundColor Cyan
Write-Host "  codon-optimize --help" -ForegroundColor Yellow
Write-Host "  python -m codon_optimizer.cli.main --help" -ForegroundColor Yellow

