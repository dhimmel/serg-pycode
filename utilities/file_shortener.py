import argparse
import os

def get_shortened_file_name(file_name, max_chars=143, verbose=False):
	"""Shorten file names to max_chars while preserving extensions."""
	if len(file_name) <= max_chars:
		return file_name
	try:
		prefix, extension = file_name.rsplit('.', 1)
		extension = '.' + extension
		if len(extension) >= max_chars:
			raise ValueError
		new_file_name = prefix[:max_chars - len(extension)] + extension
	except ValueError:
		new_file_name = file_name[:max_chars]
	if verbose:
		print file_name, 'shortened to:', new_file_name
	return new_file_name

def shorten_file_names(root_path, max_chars, verbose, no_act, force):
	"""Shorten all file names which are descendants of root_path."""
	for root, dirs, files in os.walk(root_path):
		for file_name in files:
			new_file_name = get_shortened_file_name(file_name, max_chars, verbose)
			if file_name == new_file_name:
				continue
			old_path = os.path.join(root, file_name)
			new_path = os.path.join(root, new_file_name)
			if not force and os.path.exists(new_path):
				print new_path
				raise OSError('File exists')
			if not no_act:
				os.rename(old_path, new_path)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Shorten file names')
	parser.add_argument('--max-chars', '-m', type=int, default=143)
	parser.add_argument('--verbose', '-v', action='store_true')
	parser.add_argument('--no-act', '-n', action='store_true')
	parser.add_argument('--force', '-f', action='store_true')
	parser.add_argument('--path', '-p', default=os.getcwd())
	args = parser.parse_args()
	shorten_file_names(args.path, args.max_chars, args.verbose, args.no_act, args.force)



# eCryptfs character limit: 143
# https://bugs.launchpad.net/ecryptfs/+bug/344878


