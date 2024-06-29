from setuptools import setup, find_packages

setup(
    name='GINA',
    version='0.01',
    author='Ehsan Bhatti',
    author_email='ENB360@student.bham.ac.uk',
    packages=find_packages(),
    package_data={
      'util_files': ['util_files/*'],
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'scikit-learn',
        'drugstone',
        'goatools',
        'mygene',
        'rpy2',
        'setuptools-git'
    ]

)
