from setuptools import setup

exec(open('mirror_seq/version.py').read())

INSTALL_REQUIRES = [
    'pandas==0.18.0',
    'numexpr==2.5.2',
    'pysam==0.9.0',
    'cutadapt==1.9.1',
    'tables==3.2.2',
]

setup(name='mirror_seq',
    version=__version__,
    description='The bioinformatics tool for Mirror-seq.',
    url='https://github.com/Zymo-Research/mirror-seq',
    author='Hunter Chung',
    author_email='b89603112@gmail.com',
    licence='Apache License 2.0',
    scripts=['bin/mirror-seq'],
    packages=['mirror_seq'],
    install_requires=INSTALL_REQUIRES,
    keywords='mirror sequencing next-gen hydroxymethylation bisulfite bioinformatics')
