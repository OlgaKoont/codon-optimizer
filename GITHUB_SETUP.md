# Подготовка к публикации на GitHub

## Что будет включено в репозиторий

### Основные файлы проекта:
- `codon_optimizer/` - весь исходный код
- `tests/` - тесты
- `setup.py` - установочный скрипт
- `requirements.txt` - зависимости
- `README.md` - документация

### Конфигурационные файлы:
- `.gitignore` - исключения для Git
- `.gitattributes` - настройки Git

## Что будет исключено (.gitignore)

- `DnaChisel/` - референсная реализация (не наш код)
- `__pycache__/` - кэш Python
- `*.egg-info/` - метаданные установки
- `kodon_opt/`, `venv/` - виртуальные окружения
- `*.log` - логи
- `*.exe` - установщики
- Временные файлы

## Команды для GitHub

### Инициализация репозитория (если еще не инициализирован):

```bash
git init
git add .gitignore .gitattributes
git add codon_optimizer/
git add tests/
git add setup.py
git add requirements.txt
git add README.md
git commit -m "Initial commit: Codon Optimization System"
```

### Добавление файлов (исключая DnaChisel):

```bash
# Добавить все файлы, кроме исключенных в .gitignore
git add .

# Проверить, что DnaChisel не добавлен
git status
```

### Создание репозитория на GitHub:

1. Создайте новый репозиторий на GitHub
2. Добавьте remote:
   ```bash
   git remote add origin https://github.com/yourusername/codon-optimizer.git
   ```
3. Загрузите код:
   ```bash
   git push -u origin main
   ```

## Структура репозитория

```
codon-optimizer/
├── .gitignore
├── .gitattributes
├── README.md
├── requirements.txt
├── setup.py
├── codon_optimizer/
│   ├── __init__.py
│   ├── analysis/
│   ├── antibody/
│   ├── cli/
│   ├── constraints/
│   ├── core/
│   ├── io/
│   └── specifications/
└── tests/
```

## Проверка перед коммитом

```bash
# Проверить статус
git status

# Проверить, что DnaChisel не включен
git ls-files | grep -i dnachisel

# Должно быть пусто (нет вывода)
```

