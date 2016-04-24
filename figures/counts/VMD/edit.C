for i in [filename]; do sips -s format [image type] $i --out [destination]/$i.[extension];done

for i in *.png; do sips -s format pdf $i --out $i.pdf;done

for i in *; do mv "$i" "$i.pdf"; done

for f in *.png; do mv "$f" "${f%.*}"; done