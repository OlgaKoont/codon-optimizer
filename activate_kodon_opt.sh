#!/bin/bash
# Скрипт для активации окружения kodon_opt в WSL
# Запустите: source activate_kodon_opt.sh
# Или: . activate_kodon_opt.sh

# Определяем путь к conda
CONDA_BASE="$HOME/miniconda3"
if [ ! -d "$CONDA_BASE" ]; then
    CONDA_BASE="$HOME/anaconda3"
fi

# Инициализация conda (если еще не инициализирована)
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    . "$CONDA_BASE/etc/profile.d/conda.sh"
    echo "✓ Conda инициализирована из $CONDA_BASE"
elif command -v conda &> /dev/null; then
    echo "✓ Conda уже доступна в PATH"
else
    echo "✗ Conda не найдена. Установите Miniconda:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh -b -p \$HOME/miniconda3"
    echo "  echo 'export PATH=\"\$HOME/miniconda3/bin:\$PATH\"' >> ~/.bashrc"
    return 1 2>/dev/null || exit 1
fi

# Активируем окружение
conda activate kodon_opt

if [ $? -eq 0 ]; then
    echo "✓ Окружение kodon_opt активировано!"
    echo ""
    echo "Теперь можно использовать:"
    echo "  codon-optimize --help"
    echo "  codon-optimize optimize -i input.fasta -o output.fasta --report report.html"
    echo "  python -m codon_optimizer.cli.main --help"
    echo ""
    echo "Текущая директория: $(pwd)"
else
    echo "✗ Ошибка активации окружения kodon_opt"
    echo ""
    echo "Создайте окружение:"
    echo "  conda create -n kodon_opt python=3.8 -y"
    echo "  conda activate kodon_opt"
    echo "  pip install -r requirements.txt"
    echo "  pip install -e ."
    return 1 2>/dev/null || exit 1
fi

