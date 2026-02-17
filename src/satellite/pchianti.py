import os
import tarfile
import urllib.request
from pathlib import Path


CHIANTI_URL = "https://download.chiantidatabase.org/CHIANTI_11.0.2_database.tar.gz"


def is_non_empty_dir(p: Path) -> bool:
    return p.is_dir() and any(p.iterdir())


def find_chianti_root(extract_dir: Path) -> Path:
    """
    Return the directory that should be used as XUVTOP (CHIANTI database root).
    For CHIANTI v11, the root contains (at least) 'masterlist' and 'VERSION'.
    """
    # Case A: extract_dir itself is the root
    if (extract_dir / "masterlist").exists() and (extract_dir / "VERSION").exists():
        return extract_dir

    # Case B: tarball created one top-level folder inside extract_dir
    children = [c for c in extract_dir.iterdir() if c.is_dir()]
    if len(children) == 1:
        c = children[0]
        if (c / "masterlist").exists() and (c / "VERSION").exists():
            return c

    # Case C: search one level deep (cheap & usually enough)
    for c in children:
        if (c / "masterlist").exists() and (c / "VERSION").exists():
            return c

    raise RuntimeError(
        f"Could not locate CHIANTI root under {extract_dir}. "
        f"Expected to find 'masterlist' and 'VERSION'."
    )


def prepareChianti(dir: str, logger=None) -> str:
    """
    Ensures CHIANTI DB is present and sets XUVTOP for the current Python process.

    Steps:
      1) check if dir exists and is non-empty
      2) if missing/empty, download + extract CHIANTI into dir
      3) set os.environ['XUVTOP'] to the detected CHIANTI root

    Returns:
      The path used for XUVTOP (string).
    """
    target = Path(dir).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)

    # If the directory is empty, populate it; if non-empty, we assume it might already contain CHIANTI.
    if not is_non_empty_dir(target):
        archive_path = target / Path(CHIANTI_URL).name

        # Download
        try:
            if logger:
                logger.info(f"Downloading Chianti database to {dir} from {CHIANTI_URL}")
            urllib.request.urlretrieve(CHIANTI_URL, archive_path)
        except Exception as e:
            raise RuntimeError(
                f"Failed to download CHIANTI from {CHIANTI_URL}: {e}"
            ) from e

        # Extract
        try:
            with tarfile.open(archive_path, "r:gz") as tf:
                tf.extractall(path=target)
        except Exception as e:
            raise RuntimeError(f"Failed to extract {archive_path}: {e}") from e
        finally:
            # Optional: remove archive to save space
            try:
                archive_path.unlink(missing_ok=True)
            except Exception:
                pass

    # Detect correct CHIANTI root (dir itself or a single subfolder)
    root = find_chianti_root(target)
    if logger:
        logger.info(f"Setting Chianti db path to {root}")

    # "Export" for this Python process (and anything it spawns)
    os.environ["XUVTOP"] = str(root)

    # Optional sanity check (won't error if you don't have pyneb installed yet)
    # You can uncomment if you want it to fail fast when ChiantiPy/PyNeb can't see it.
    #
    # import pyneb as pn
    # _ = pn.utils.pn_chianti.Chianti_getA("o_3")

    return os.environ["XUVTOP"]
