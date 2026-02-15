# Codon Optimization System

Система кодон-оптимизации для производства антител в клетках CHO с мультикритериальной оптимизацией.

## Установка окружения

### Требования
- Python 3.8 или выше
- pip

### Установка

1. **Создайте виртуальное окружение:**

```bash
# Windows
python -m venv kodon_opt
kodon_opt\Scripts\activate

# Linux/Mac
python -m venv kodon_opt
source kodon_opt/bin/activate
```

2. **Установите зависимости:**

```bash
pip install -r requirements.txt
pip install -e .
```

3. **Проверьте установку:**

```bash
codon-optimize --help
```

## Запуск приложения

### Важно
Система принимает **белковые (аминокислотные) последовательности** на входе. Если подана ДНК-последовательность, она автоматически транслируется в белок.

### Формат входного файла

Входной FASTA файл должен содержать белковую последовательность:

```
>protein_sequence
EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISSGGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR...
```

### 1. Базовая оптимизация (одна последовательность)

```bash
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --host CHO \
    --report report.html
```

### 2. Генерация нескольких решений (по разным критериям)

Генерирует до 5 решений, каждое оптимизировано по разным критериям:
- **balanced** - лучший общий score (сбалансированный)
- **high_cai** - максимальный CAI (codon adaptation index)
- **optimal_gc** - оптимальное содержание GC
- **minimal_motifs** - минимум проблемных мотивов
- **best_codon_pairs** - лучшее использование пар кодонов

```bash
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --num-solutions 5 \
    --host CHO \
    --report report.html
```

### 3. Pareto-оптимизация (недоминируемые решения)

```bash
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --mode pareto \
    --pareto-solutions 10 \
    --host CHO \
    --report report.html
```

### 4. Оптимизация антител (тяжелая и легкая цепи)

```bash
codon-optimize optimize \
    --input input.fasta \
    --output optimized.fasta \
    --mode antibody \
    --heavy-chain HC.fasta \
    --light-chain LC.fasta \
    --host CHO \
    --report report.html
```

### 5. Оптимизация с ограничениями

```bash
# Удаление рестриктазных сайтов
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --restriction-sites "BsaI_site,BsmBI_site,NotI_site" \
    --host CHO \
    --report report.html

# Целевой диапазон экспрессии
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --target-expression 0.7-0.9 \
    --host CHO \
    --report report.html
```

### 6. Анализ без оптимизации

```bash
codon-optimize analyze \
    --input protein.fasta \
    --output analysis.html
```

### Параметры командной строки

| Параметр | Описание | По умолчанию |
|----------|----------|--------------|
| `--input`, `-i` | Входной FASTA файл (белок или ДНК) | Обязательно |
| `--output`, `-o` | Выходной FASTA файл | Обязательно |
| `--host` | Организм-хозяин | CHO |
| `--mode` | Режим: `standard`, `pareto`, `antibody` | standard |
| `--num-solutions` | Количество решений (1-5) | 1 |
| `--pareto-solutions` | Количество Pareto решений | 10 |
| `--restriction-sites` | Список рестриктазных сайтов для удаления | - |
| `--target-expression` | Целевой диапазон экспрессии (0.7-0.9) | - |
| `--report` | Путь к отчету | - |
| `--report-format` | Формат отчета: `text`, `json`, `html` | html |

## Чтение результатов

### 1. Выходной FASTA файл

#### Одна последовательность
```
>protein_sequence_optimized
ATGGAA...
```

#### Несколько последовательностей
```
>optimized_balanced
ATGGAA...
>optimized_high_cai
ATGGAG...
>optimized_optimal_gc
ATGGAC...
>optimized_minimal_motifs
ATGGAT...
>optimized_best_codon_pairs
ATGGAC...
```

### 2. HTML отчет

Откройте файл отчета в браузере. Отчет содержит:

#### Для одной последовательности:
- **Общая информация**: длина последовательности, метод оптимизации
- **Scores**: общий score, CAI, GC, пары кодонов, структура мРНК, мотивы
- **CAI (Codon Adaptation Index)**: значение CAI и распределение кодонов
- **GC Content**: общее содержание GC, профиль по последовательности
- **mRNA Structure**: структура 5'-конца, стабильные структуры
- **Motifs**: найденные проблемные мотивы (TATA-box, cryptic promoters и т.д.)
- **Codon Usage**: использование кодонов, избежание истощения тРНК
- **CQA Prediction** (для антител): предсказание гликозилирования, агрегации, зарядовых вариантов

#### Для нескольких последовательностей:
- **Сравнительная таблица** всех решений
- **Детальный анализ** каждого решения
- **Рекомендации** по выбору решения

### 3. JSON отчет

Структурированные данные для программной обработки:

```json
{
  "sequences": [
    {
      "sequence_id": "optimized_balanced",
      "dna_sequence": "ATGGAA...",
      "scores": {
        "total_score": 0.9041,
        "cai_score": 0.8234,
        "gc_score": 0.9123,
        ...
      },
      "metadata": {
        "strategy": "balanced",
        "cai": 0.8234,
        "gc_content": 0.5123,
        ...
      }
    }
  ]
}
```

### 4. Текстовый отчет

Человекочитаемый текстовый формат с детальной информацией по каждому решению.

## Интерпретация результатов

### Общий Score (0.0 - 1.0)
- **> 0.8**: Отличная оптимизация
- **0.6 - 0.8**: Хорошая оптимизация
- **< 0.6**: Требует улучшения

### CAI (Codon Adaptation Index)
- **> 0.8**: Высокий CAI, хорошая адаптация к хосту
- **0.6 - 0.8**: Средний CAI
- **< 0.6**: Низкий CAI, возможны проблемы с экспрессией

### GC Content
- **0.3 - 0.8**: Нормальный диапазон для CHO
- **< 0.3 или > 0.8**: Могут быть проблемы со стабильностью мРНК

### Motifs
- **0 мотивов**: Идеально
- **1-5 мотивов**: Приемлемо
- **> 5 мотивов**: Могут быть регуляторные проблемы

### Выбор решения

При генерации нескольких решений:

1. **Для максимальной экспрессии**: выберите `optimized_high_cai`
2. **Для избежания проблем**: выберите `optimized_minimal_motifs`
3. **Для общего использования**: выберите `optimized_balanced`
4. **Для оптимального GC**: выберите `optimized_optimal_gc`

## Примеры использования

### Пример 1: Базовая оптимизация

```bash
codon-optimize optimize \
    --input IgG1_HC_example.fasta \
    --output IgG1_HC_optimized.fasta \
    --host CHO \
    --report IgG1_HC_report.html
```

### Пример 2: Несколько решений

```bash
codon-optimize optimize \
    --input IgG1_HC_example.fasta \
    --output IgG1_HC_multiple.fasta \
    --num-solutions 5 \
    --host CHO \
    --report IgG1_HC_multiple_report.html
```

### Пример 3: С ограничениями

```bash
codon-optimize optimize \
    --input protein.fasta \
    --output optimized.fasta \
    --restriction-sites "BsaI_site,BsmBI_site" \
    --target-expression 0.75-0.85 \
    --host CHO \
    --report report.html
```

## Структура проекта

```
codon_optimizer/
├── core/              # Ядро оптимизации
├── analysis/          # Анализ последовательностей
├── constraints/       # Ограничения
├── specifications/    # Спецификации (constraints/objectives)
├── antibody/         # Специфичные для антител функции
├── io/               # Ввод/вывод и отчеты
└── cli/              # Командная строка
```

## Поддержка

При возникновении проблем проверьте:
1. Активировано ли виртуальное окружение
2. Установлены ли все зависимости (`pip install -r requirements.txt`)
3. Правильный ли формат входного файла (белковая последовательность)
4. Логи в консоли для диагностики ошибок
