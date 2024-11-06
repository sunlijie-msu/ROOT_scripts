from gtts import gTTS

# Replace 'Your text here' with the text you want to convert
text = '''Your text here'''

# Set the language (e.g., 'en' for English)
language = 'en'

# Create the gTTS object
speech = gTTS(text=text, lang=language, slow=False)

# Save the converted audio to a file
speech.save("output.mp3")

print("Audio content written to output.mp3")
