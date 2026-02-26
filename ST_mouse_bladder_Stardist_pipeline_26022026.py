
pip install stardist csbdeep tifffile scikit-image matplotlib pandas numpy opencv-python

import os
import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt

from stardist.models import StarDist2D
from csbdeep.utils import normalize
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries


# ============================
# USER PARAMETERS
# ============================

image_path = "input_image.tif"
transcript_file = "transcripts.csv"
output_dir = "stardist_output"

model_name = "2D_versatile_fluo"

os.makedirs(output_dir, exist_ok=True)


# ============================
# Load image
# ============================

img = tifffile.imread(image_path)
img_norm = normalize(img, 1, 99.8)


# ============================
# Run StarDist
# ============================

model = StarDist2D.from_pretrained(model_name)
labels, details = model.predict_instances(img_norm)


# ============================
# Save segmentation outputs
# ============================

tifffile.imwrite(
    os.path.join(output_dir, "nuclei_labels.tif"),
    labels.astype(np.uint16)
)

boundaries = find_boundaries(labels)
overlay = np.dstack([img_norm]*3)
overlay[boundaries] = [1, 0, 0]

plt.imshow(overlay)
plt.axis("off")
plt.savefig(
    os.path.join(output_dir, "nuclei_overlay.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()


# ============================
# Extract centroids
# ============================

regions = regionprops(labels)

coords = np.array([r.centroid for r in regions])
coords_df = pd.DataFrame(coords, columns=["y", "x"])
coords_df.index.name = "nucleus_id"

coords_df.to_csv(
    os.path.join(output_dir, "banksy_coords.csv")
)


# ============================
# Assign transcripts to nuclei
# ============================

transcripts = pd.read_csv(transcript_file)

counts = {}

for region in regions:
    label_id = region.label
    mask = labels == label_id

    # assume transcripts have integer pixel coords
    hits = transcripts[
        mask[
            transcripts["y"].astype(int),
            transcripts["x"].astype(int)
        ]
    ]

    gene_counts = hits["gene"].value_counts()
    counts[label_id] = gene_counts

counts_df = pd.DataFrame(counts).fillna(0).astype(int)
counts_df.index.name = "gene"

counts_df.to_csv(
    os.path.join(output_dir, "banksy_counts.csv")
)

print("Done.")
