import cv2
import os
import re

def natural_sort_key(filename):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', filename)]

image_folder = 'Results/Geometry_study'
output_video = 'output.mp4'

images = sorted([img for img in os.listdir(image_folder) if img.endswith(('.png', '.jpg', '.jpeg'))],key=natural_sort_key)

first_image_path = os.path.join(image_folder, images[0])
first_image = cv2.imread(first_image_path)
height, width, _ = first_image.shape

fps = 30
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

for image in images:
    img_path = os.path.join(image_folder, image)
    frame = cv2.imread(img_path)
    
    video.write(frame)

video.release()
cv2.destroyAllWindows()

print("Video saved as", output_video)


