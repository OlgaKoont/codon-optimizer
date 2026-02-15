#!/bin/bash
# Скрипт для полной настройки окружения в WSL
# Запустите: bash setup_wsl.sh

set -e  # Остановка при ошибке

echo "=========================================="
echo "Настройка codon-optimizer в WSL"
echo "=========================================="
echo ""

# Определяем путь к проекту
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Директория проекта: $PROJECT_DIR"
cd "$PROJECT_DIR"

# Проверяем наличие conda
if ! command -v conda &> /dev/null; then
    echo ""
    echo "Conda не найдена. Установка Miniconda..."
    echo ""
    
    cd ~
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
    
    # Добавляем conda в PATH
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    
    # Инициализируем conda
    $HOME/miniconda3/bin/conda init bash
    
    echo "✓ Miniconda установлена"
    echo "  Перезапустите терминал или выполните: source ~/.bashrc"
    echo ""
    
    # Используем conda из установленного пути
    export PATH="$HOME/miniconda3/bin:$PATH"
    source $HOME/miniconda3/etc/profile.d/conda.sh
else
    echo "✓ Conda найдена: $(which conda)"
    # Инициализируем conda если нужно
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source $HOME/miniconda3/etc/profile.d/conda.sh
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source $HOME/anaconda3/etc/profile.d/conda.sh
    fi
fi

echo ""
echo "Создание conda окружения kodon_opt..."
echo ""

# Создаем окружение если его нет
if conda env list | grep -q "kodon_opt"; then
    echo "✓ Окружение kodon_opt уже существует"
else
    conda create -n kodon_opt python=3.8 -y
    echo "✓ Окружение kodon_opt создано"
fi

# Активируем окружение
conda activate kodon_opt

echo ""
echo "Установка зависимостей..."
echo ""

# Устанавливаем зависимости
pip install --upgrade pip
pip install -r requirements.txt
pip install -e .

echo ""
echo "Проверка установки..."
echo ""

# Проверяем установку
python --version
echo ""
echo "Установленные пакеты:"
pip list | grep -E "biopython|numpy|pymoo|click|codon-optimizer" || true

echo ""
echo "=========================================="
echo "✓ Установка завершена!"
echo "=========================================="
echo ""
echo "Для активации окружения выполните:"
echo "  source activate_kodon_opt.sh"
echo ""
echo "Или вручную:"
echo "  conda activate kodon_opt"
echo ""
echo "Запуск оптимизации:"
echo "  codon-optimize optimize -i IgG1_HC_example.fasta -o output.fasta --report report.html"
echo ""

