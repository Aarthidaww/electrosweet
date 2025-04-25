import multiprocessing
import signal
import sys
from app import app  # Replace with your actual app import

# Configuration
workers = multiprocessing.cpu_count() * 2 + 1
worker_class = "gthread"
threads = 4
bind = "0.0.0.0:10000"
timeout = 120
keepalive = 5
graceful_timeout = 30

def on_exit(server):
    # Add custom cleanup logic here
    print("Performing pre-shutdown cleanup...")

def worker_int(worker):
    print(f"Worker {worker.pid} received interrupt")
    worker.cancel_init()

def worker_abort(worker):
    print(f"Worker {worker.pid} received abort signal")

# Handle OS signals
def handle_term(signum, frame):
    print("SIGTERM received: Starting graceful shutdown")
    sys.exit(0)

signal.signal(signal.SIGTERM, handle_term)
