import logging
import sys

class PyNebLogRedirector:
    """Redirects only PyNeb's stderr output to the logger while keeping the rest in stderr."""
    def __init__(self, logger, level=logging.WARNING):
        self.logger = logger
        self.level = level
        self.stderr = sys.stderr  # Store original stderr

    def write(self, message):
        if ("pyneb" in message.lower()) or (message.lower().startswith("warng")):  # Only capture PyNeb messages
            self.logger.log(self.level, message.strip())
        else:
            self.stderr.write(message)  # Pass non-PyNeb messages to original stderr

    def flush(self):  # Needed for compatibility with sys.stderr
        self.stderr.flush()

def setup_logger(name: str, log_type=logging.INFO, log_to_file=None):
    """
    Parameter:
        name (str): name of logger
        log_type (logging.TYPE): TYPE can be any of:
            logging.DEBUG → Logs everything (debug, info, warning, etc.).
            logging.INFO → Logs info, warning, error, critical (but not debug).
            logging.WARNING → Logs only warnings, errors, and critical messages.
            logging.ERROR → Logs only errors and critical messages.
            logging.CRITICAL → Logs only critical messages.
        log_to_file (str): If not None, then the logger will direct all 
          messages to a file named log_to_file (which will be created at
          first function call).
    """
    logger = logging.getLogger(name)
    logger.setLevel(log_type)
    # Define log format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    # Create a log handler (console or file based on user choice)
    if log_to_file is not None:
        handler = logging.FileHandler(log_to_file)
    else:
        handler = logging.StreamHandler()
    handler.setLevel(log_type)
    handler.setFormatter(formatter)
    # Add handler to logger if not already added (prevent duplicates)
    if not logger.hasHandlers():
        logger.addHandler(handler)
    return logger
