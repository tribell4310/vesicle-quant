"""

preprocess.py

This is built from VYA's jupyter notebook.  Goal is to quantify membrane area in EM images.

"""


# Import necessary packages: It is best if you have the packages installed in a conda environment.
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("Agg")
mpl.rc("figure", max_open_warning = 0)
from scipy.ndimage import gaussian_filter
import skimage
from skimage import data
from skimage import img_as_float
from skimage.feature import peak_local_max
from skimage.morphology import reconstruction
from skimage.color import rgb2gray
import cv2 as cv
from scipy import ndimage as ndi
from skimage import io
from datetime import datetime


def main():
	start_time = datetime.now()
	# Load the images from local directory:
	inFiles = []
	for item in os.listdir(os.getcwd()):
		if item.endswith(".tif") or item.endswith(".TIF"):
			inFiles.append(item)

	with open("membrane_counts.csv", "w") as g:
		g.write("FILE,MEMBRANE PIXELS,MEMBRANE AREA AS FRACTION OF IMAGE\n")
		print("\nProcessing tif images, this could take a while.  Go make a cup of tea.\n")
		for i, each_file in enumerate(inFiles):
			print("\t"+str(i+1)+" / "+str(len(inFiles))+" ("+each_file+") ... ")

			image1 = img_as_float(io.imread(each_file))
			image1_gaus = gaussian_filter(image1, 1)

			seed = np.copy(image1)
			seed[1:-1, 1:-1] = image1.min()
			mask1 = image1

			dilated1 = reconstruction(seed, mask1, method='dilation')

			# Subtracting the dilated image leaves an image with just the ROIs and a flat, black background, as shown below.
			final_image1 = image1 - dilated1

			# Otsu thresholding
			otsu_thresh = otsu_threshhold(final_image1)

			# Denoising step via distance transform
			no_salt = desalt_binary_mask(otsu_thresh.copy(), 2)

			# Total pixel count
			total_px = final_image1.shape[0] * final_image1.shape[1]
	
			# Apply threshold
			pixel_count = count_pixels(no_salt)
			g.write(str(each_file)+","+str(pixel_count)+","+str(pixel_count / total_px)+"\n")

			# Final plot out
			fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(nrows=1, ncols=5, figsize=(24, 4), sharex=True, sharey=True)

			ax0.imshow(image1, cmap='gray')
			ax0.set_title('Raw Image')
			ax0.axis('off')

			ax1.imshow(dilated1, vmin=image1.min(), vmax=image1.max(), cmap='gray')
			ax1.set_title('Dilated')
			ax1.axis('off')

			ax2.imshow(final_image1, cmap='gray')
			ax2.set_title('image1 - dilated1')
			ax2.axis('off')

			ax3.imshow(otsu_thresh, cmap="gray")
			ax3.set_title("otsu")
			ax3.axis('off')

			ax4.imshow(no_salt, cmap="gray")
			ax4.set_title("desalter")
			ax4.axis('off')

			fig.tight_layout()
			plt.savefig(no_ext(each_file)+"_processing.png")
			plt.clf()

	# Output
	end_time = datetime.now()
	delta = end_time - start_time
	print("\nTotal images processed:\t"+str(len(inFiles)))
	print("Total run time:\t\t"+str(delta))
	print("Time per micrograph:\t"+str(delta / len(inFiles)))


def desalt_binary_mask(in_img, dt_bound):
	pass
	# Get a distance transform of the binary mask
	dist_transform = ndi.distance_transform_edt(in_img)

	# Heuristic: if the dist_transform pixel is <= 3 AND no neighbor is <= 3, then flip to black.
	for i in range(0, in_img.shape[0]):
		for j in range(0, in_img.shape[1]):
			if dist_transform[i][j] <= dt_bound:
				if check_eight_neighbors(dist_transform, i, j, in_img.shape, dt_bound) == False:
					in_img = flip_px_to_zero(in_img, i, j)

	return in_img



def check_eight_neighbors(dt, x, y, shape, bound):
	""" This will break if the image is non-dimensional, but then again it should have broken long before that then. """

	# This is a nine-part switch statement softcoded for all possible pixel values.
	# Returns False if all neighbor pixels are <= bound value, True if any pixel exceeds the bound value.

	shape_x = shape[0]
	shape_y = shape[1]

	if x == 0 and y == 0: # corner - 3 neighbors
		if (dt[x][y+1] <= bound) and (dt[x+1][y] <= bound) and (dt[x+1][y+1] <= bound):
			return False
		else:
			return True

	elif x == 0 and y == shape_y-1: # corner - 3 neighbors
		if (dt[x][y-1] <= bound) and (dt[x+1][y-1] <= bound) and (dt[x+1][y] <= bound):
			return False
		else:
			return True

	elif x == shape_x-1 and y == 0: # corner - 3 neighbors
		if (dt[x-1][y] <= bound) and (dt[x-1][y+1] <= bound) and (dt[x][y+1] <= bound):
			return False
		else:
			return True

	elif x == shape_x-1 and y == shape_y-1: # corner - 3 neighbors
		if (dt[x-1][y-1] <= bound) and (dt[x-1][y] <= bound) and (dt[x][y-1] <= bound):
			return False
		else:
			return True

	elif x == 0: # edge - 5 neighbors
		if (dt[x][y-1] <= bound) and (dt[x][y+1] <= bound) and (dt[x+1][y-1] <= bound) and (dt[x+1][y] <= bound) and (dt[x+1][y+1] <= bound):
			return False
		else:
			return True

	elif x == shape_x-1: # edge - 5 neighbors
		if (dt[x-1][y-1] <= bound) and (dt[x-1][y] <= bound) and (dt[x-1][y+1] <= bound) and (dt[x][y-1] <= bound) and (dt[x][y+1] <= bound):
			return False
		else:
			return True

	elif y == 0: # edge - 5 neighbors
		if (dt[x-1][y] <= bound) and (dt[x+1][y] <= bound) and (dt[x-1][y+1] <= bound) and (dt[x][y+1] <= bound) and (dt[x+1][y+1] <= bound):
			return False
		else:
			return True

	elif y == shape_y-1: # edge - 5 neighbors
		if (dt[x-1][y-1] <= bound) and (dt[x][y-1] <= bound) and (dt[x+1][y-1] <= bound) and (dt[x-1][y] <= bound) and (dt[x+1][y] <= bound):
			return False
		else:
			return True

	else: # internal - 8 neighbors
		if (dt[x-1][y-1] <= bound) and (dt[x-1][y] <= bound) and (dt[x-1][y+1] <= bound) and (dt[x][y-1] <= bound) and (dt[x][y+1] <= bound) and (dt[x+1][y-1] <= bound) and (dt[x+1][y] <= bound) and (dt[x+1][y+1] <= bound):
			return False
		else:
			return True



def simple_watershed(img):
	# Calculate distance transform of binary mask and find local maxima
	print("segmentin")
	dist_transform = ndi.distance_transform_edt(img)
	local_max_identities = peak_local_max(dist_transform, min_distance=1)
	local_max_boolean = np.zeros((img.shape[0], img.shape[1]), dtype="uint16")
	for local_max in local_max_identities:
		local_max_boolean[local_max[0]][local_max[1]] = 1

	# Segment each local maximum into a marked object and extract its properties
	markers, _ = ndi.label(local_max_boolean)
	segmented = skimage.segmentation.watershed(255-dist_transform, markers, mask=img, compactness=0.005)
	object_labels = skimage.measure.label(segmented)
	properties = skimage.measure.regionprops(object_labels)
	print(properties)
	exit()

	return segmented


def background_mask_dt(inMask):
	# Define some static variables to avoid unnecessary calls
	mask_shape_x = inMask.shape[0]
	mask_shape_y = inMask.shape[1]

	# Create the new mask object
	bg_mask = np.copy(inMask) # white circles on black bg

	# Calculate a distance transform of bg_mask
	bg_mask_dt = ndi.distance_transform_edt(bg_mask)

	return bg_mask


def flip_px_to_zero(mask, x, y):
	mask[x, y] = 0
	return mask


def otsu_threshhold(inImg):
	# Actually perform threshholding.
	inImg = inImg * 255 / np.max(inImg)
	otsu_threshold, image_result = cv.threshold(inImg.astype("uint8"), 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
	
	return image_result


def count_pixels(in_mask):
	counter = 0
	for i in range(0, in_mask.shape[0]):
		for j in range(0, in_mask.shape[1]):
			if in_mask[i][j] > 0:
				counter += 1

	return counter


def janky_heuristic_mask(in_img, thresh):
	px_count = 0
	out_img = np.zeros((in_img.shape[0], in_img.shape[1]))
	for i in range(0, in_img.shape[0]):
		for j in range(0, in_img.shape[1]):
			if in_img[i][j] > thresh:
				out_img[i][j] = in_img[i][j]
				px_count += 1
	return out_img, px_count


def clean_large_numbers(inVal):
	inStr = str(inVal)
	if len(inStr) <= 3:
		return inStr
	else:
		outStr = ""
		while len(inStr) > 3:
			outStr = inStr[-3:] + "," + outStr
			inStr = inStr[:-3]
		if len(inStr) != 0:
			outStr = inStr + "," + outStr
		return outStr[:-1]


def no_ext(inStr):
	"""
	Takes an input filename and returns a string with the file extension removed.

	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find(".", prevPos+1)
	return inStr[0:prevPos]


if __name__ == "__main__":
	if len(sys.argv) == 1:
		main()
	else:
		print("Check inputs: foo.py")