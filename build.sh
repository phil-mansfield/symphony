rm -r build/* dist/*
python3 setup.py sdist bdist_wheel
python3 -m twine upload --repository symlib dist/*
