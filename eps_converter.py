from PIL import Image

# Open the image file
img = Image.open("D:\\S\\Picture\\Experimental photos\\Photoshop Detector\\Photoshop NPR\\PXCT_cartoon.jpg")

# Convert image to RGB if it's not already
if img.mode != 'RGB':
    img = img.convert('RGB')

# Set the DPI (try a value like 150 or 300)
dpi = (20, 20)

# Save the image as EPS with the specified DPI
img.save("D:\\S\\Picture\\Experimental photos\\Photoshop Detector\\Photoshop NPR\\PXCT_cartoon.eps", format='EPS', dpi=dpi)

print("Conversion complete.")
