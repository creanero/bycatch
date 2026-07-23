"""Per-target configuration template for HAT-P-32b.

Usage pattern:
1) Fill TARGET_POSITIONS once you have coordinates.
2) In scripts, import from this profile instead of config.py.
"""

from config import *  # noqa: F401,F403

# Override only target-specific values.
RAW_DATA_FOLDER = "raw_images_hatp32"
ALIGNED_FOLDER = "aligned_images_hatp32"

# TODO: Replace these placeholders with real (x, y) tuples.
# Convention: target first, then comparison stars.
TARGET_POSITIONS = (
    (0.0, 0.0),  # TODO target star (HAT-P-32b)
    (0.0, 0.0),  # TODO comparison star 1
    (0.0, 0.0),  # TODO comparison star 2
)
