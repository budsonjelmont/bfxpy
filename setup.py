import setuptools

setuptools.setup(
    name='bfxpy',
    version='1.0.0',
    url='https://github.com/budsonjelmont/bfxpy',
    author='Judson Belmont',
    author_email='budsonjelmont@gmail.com',
    description='Various bioinformatics utilities',
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=['pandas','numpy','pyvcf'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ]
)

