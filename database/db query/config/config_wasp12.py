"""Per-target configuration template for WASP-12-b.

Usage pattern:
1) Fill TARGET_POSITIONS once you have coordinates.
2) In scripts, import from this profile instead of config.py.
"""

from config import *  # noqa: F401,F403

# Override only target-specific values.
RAW_DATA_FOLDER = "raw_images_wasp122"
ALIGNED_FOLDER = "aligned_images_wasp12"

# TODO: Replace these placeholders with real (x, y) tuples.
# Convention: target first, then comparison stars.
TARGET_POSITIONS = (
    (0.0, 0.0),  # TODO target star (WASP-12-b)
    (0.0, 0.0),  # TODO comparison star 1
    (0.0, 0.0),  # TODO comparison star 2
)
