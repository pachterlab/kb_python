from setuptools import setup

setup(
    name='kb_python',
    version='0.0.0',
    description='',
    url='https://github.com/pachterlab/kb_python',
    author='Kyung Hoi (Joseph) Min',
    author_email='phoenixter96@gmail.com',
    license='MIT',
    packages=['kb_python'],
    zip_safe=False,
    entry_points={
        'console_scripts': ['kb=kb_python.main:main'],
    },
)
