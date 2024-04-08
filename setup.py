from setuptools import setup, find_packages

setup(
    name='GINA',
    version='0.01',
    url='',
    author='Ehsan Bhatti',
    author_email='',
    description='',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'scikit-learn',
        'drugstone',
        'goatools',
        'mygene',
        'rpy2'
    ]

)
