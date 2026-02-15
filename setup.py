"""
Setup script for codon optimization package.

Рекомендуемое окружение: kodon_opt

Для создания окружения выполните:
    conda create -n kodon_opt python=3.8
    conda activate kodon_opt
    pip install -r requirements.txt
    pip install -e .
"""

from setuptools import setup, find_packages

setup(
    name='codon-optimizer',
    version='0.1.0',
    description='Codon Optimization System for Antibody Production',
    author='Bioinformatics Team',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.81',
        'numpy>=1.24.0',
        'scipy>=1.10.0',
        'pandas>=2.0.0',
        'matplotlib>=3.7.0',
        'pymoo>=0.6.0',
        'scikit-learn>=1.3.0',
        'click>=8.1.0',
    ],
    entry_points={
        'console_scripts': [
            'codon-optimize=codon_optimizer.cli.main:cli',
        ],
    },
    python_requires='>=3.8',
    long_description="""
    Codon Optimization System for Antibody Production
    
    Рекомендуемое окружение: kodon_opt
    
    Установка:
    1. Создайте conda окружение:
       conda create -n kodon_opt python=3.8
       conda activate kodon_opt
    
    2. Установите зависимости:
       pip install -r requirements.txt
       pip install -e .
    
    3. Используйте систему:
       codon-optimize optimize --input sequence.fasta --output optimized.fasta
    """,
)


