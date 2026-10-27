#!/usr/bin/env python3
"""Remove local file-name keys from a compiled PDF without changing its pages."""
from pathlib import Path
import argparse, json, os, tempfile
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DictionaryObject, NameObject

parser = argparse.ArgumentParser()
parser.add_argument('pdf', type=Path)
args = parser.parse_args()
reader = PdfReader(args.pdf)
writer = PdfWriter()
writer.clone_document_from_reader(reader)
removed = 0
for obj in writer._objects:
    if isinstance(obj, DictionaryObject) and '/PTEX.FileName' in obj:
        del obj[NameObject('/PTEX.FileName')]
        removed += 1
if removed:
    fd, tmp = tempfile.mkstemp(prefix='.pdf-metadata-', suffix='.pdf', dir=args.pdf.parent)
    try:
        with os.fdopen(fd, 'wb') as output:
            writer.write(output)
        os.replace(tmp, args.pdf)
    finally:
        if Path(tmp).exists():
            Path(tmp).unlink()
print(json.dumps({'removed_local_file_name_keys': removed, 'pages': len(reader.pages)}))
