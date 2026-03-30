from setuptools import setup

setup(
    name="onyx",
    version="0.0.1",
    py_modules=["onyx"],
    install_requires=["tqdm"],
    entry_points={
        "console_scripts": [
            "onyx=onyx:main"
        ]
    }
)

