import sys
import os


assert len(sys.argv) == 4, "Bad number of arguments."
file_path = sys.argv[1]
prefix = sys.argv[2]
comment = sys.argv[3]

assert os.path.isfile(file_path)

with open(file_path) as read_handler:
	text = read_handler.read()
	sys.stdout.write(text.format(prefix=prefix, comment=comment))

