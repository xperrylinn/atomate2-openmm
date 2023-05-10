import logging

# Define logger name
logger_name = "atomate2-openmm"

# Create the logger object
logger = logging.getLogger(logger_name)
logger.setLevel(logging.DEBUG)

# Create the console handler and set its log level
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Create the file handler and set its log level
file_handler = logging.FileHandler(f"{logger_name}.log")
file_handler.setLevel(logging.DEBUG)

# Create the log format for both handlers
log_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(log_format)
file_handler.setFormatter(log_format)

# Add both handlers to the logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)
