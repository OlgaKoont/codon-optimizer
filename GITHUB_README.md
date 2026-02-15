# Инструкция по загрузке на GitHub

## Быстрый старт

### 1. Инициализация Git репозитория

```bash
cd C:\Users\466259\Desktop\Case
git init
```

### 2. Добавление файлов

```bash
# Добавить все файлы (DnaChisel автоматически исключен через .gitignore)
git add .

# Проверить статус
git status
```

### 3. Проверка, что DnaChisel не включен

```bash
# Должно быть пусто (нет вывода)
git ls-files | findstr /i dnachisel
```

### 4. Первый коммит

```bash
git commit -m "Initial commit: Codon Optimization System for Antibody Production"
```

### 5. Создание репозитория на GitHub

1. Зайдите на https://github.com
2. Создайте новый репозиторий (например, `codon-optimizer`)
3. **НЕ** инициализируйте с README, .gitignore или лицензией

### 6. Подключение и загрузка

```bash
# Добавить remote (замените YOUR_USERNAME на ваш GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/codon-optimizer.git

# Загрузить код
git branch -M main
git push -u origin main
```

## Что будет включено

✅ **Включено:**
- `codon_optimizer/` - весь исходный код
- `tests/` - тесты
- `setup.py` - установочный скрипт
- `requirements.txt` - зависимости
- `README.md` - документация
- `.gitignore` - исключения
- `.gitattributes` - настройки Git

❌ **Исключено (через .gitignore):**
- `DnaChisel/` - референсная реализация
- `__pycache__/` - кэш Python
- `*.egg-info/` - метаданные
- `kodon_opt/`, `venv/` - виртуальные окружения
- `*.log` - логи
- `*.exe` - установщики
- Временные файлы

## Структура репозитория на GitHub

```
codon-optimizer/
├── .gitignore
├── .gitattributes
├── README.md
├── requirements.txt
├── setup.py
├── codon_optimizer/
│   ├── __init__.py
│   ├── analysis/          # Анализ последовательностей
│   ├── antibody/           # Специфичные для антител функции
│   ├── cli/                # Командная строка
│   ├── constraints/        # Ограничения
│   ├── core/               # Ядро оптимизации
│   ├── io/                 # Ввод/вывод и отчеты
│   └── specifications/     # Спецификации (constraints/objectives)
└── tests/                  # Тесты
```

## Проверка перед загрузкой

```bash
# 1. Проверить статус
git status

# 2. Убедиться, что DnaChisel не включен
git ls-files | findstr /i dnachisel
# (должно быть пусто)

# 3. Посмотреть список файлов, которые будут загружены
git ls-files

# 4. Проверить размер репозитория
git count-objects -vH
```

## Дополнительные файлы (опционально)

Если хотите включить примеры результатов, раскомментируйте в `.gitignore`:
```gitignore
# *.fasta
# *.html
```

Если хотите включить скрипты активации, раскомментируйте:
```gitignore
# activate_kodon_opt.ps1
# activate_kodon_opt.sh
```

