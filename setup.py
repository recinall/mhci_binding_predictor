from setuptools import setup, find_packages

setup(
    name="mhci_binding_predictor",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'requests',
        'click',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'mhci_predict=mhci_binding_predictor.cli:main',
        ],
    },
    author="Your Name",
    description="MHC-I binding prediction and immunogenicity analysis tool",
    license="MIT"
)
