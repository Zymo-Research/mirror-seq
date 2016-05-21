from setuptools import setup

exec(open('mirror_seq/version.py').read())
LONG_DESCRIPTION = '''
Please visit the GitHub repo (https://github.com/Zymo-Research/mirror-seq) for detail information.
'''


INSTALL_REQUIRES = [
    'pandas>=0.18.0',
    'pysam>=0.9.0',
    'cutadapt==1.9.1',
]

setup(name='mirror_seq',
    version=__version__,
    description='The bioinformatics tool for Mirror-seq.',
    long_description=LONG_DESCRIPTION,
    url='https://github.com/Zymo-Research/mirror-seq',
    author='Hunter Chung',
    author_email='b89603112@gmail.com',
    licence='Apache License 2.0',
    scripts=['bin/mirror-seq', 'bin/mirror-trim', 'bin/mirror-call'],
    packages=['mirror_seq'],
    install_requires=INSTALL_REQUIRES,
    classifiers=['Programming Language :: Python :: 2.7'],
    keywords='mirror sequencing next-gen hydroxymethylation bisulfite bioinformatics')
