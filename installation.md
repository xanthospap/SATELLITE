### Clone the repository
```
git clone https://github.com/xanthospap/SATELLITE.git
```
and get into the root directory, i.e. `cd SATELLITE`.

### Fetch branches
```
git fetch --all
```

### Create a new branch to track remote/xanthos

If the branch `xanthos` does not exist:
```
git checkout -b xanthos origin/xanthos
```

### Pull changes (latest)
```
git pull origin xanthos
```

### Install via pip, using local folder
```
pip install -e .
```
