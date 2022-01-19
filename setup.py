from setuptools import find_packages, setup


def read(path):
    with open(path, 'r') as f:
        return f.read()


long_description = read('README.md')

setup(
    name='kb_python',
    version='0.27.0',
    url='https://github.com/pachterlab/kb_python',
    author='Kyung Hoi (Joseph) Min',
    author_email='phoenixter96@gmail.com',
    maintainer='Pachter Lab',
    maintainer_email='lpachter@caltech.edu',
    description='Python wrapper around kallisto | bustools for scRNA-seq analysis',  # noqa
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='kallisto bustools',
    python_requires='>=3.6',
    license='BSD',
    packages=find_packages(exclude=('tests', 'tests.*', 'docs')),
    zip_safe=False,
    include_package_data=True,
    install_requires=read('requirements.txt').strip().split('\n'),
    entry_points={
        'console_scripts': ['kb=kb_python.main:main'],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Utilities',
    ],
)
