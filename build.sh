rm build/*
python3 setup.py sdist bdist_wheel
twine upload -r testpypi dist/*
