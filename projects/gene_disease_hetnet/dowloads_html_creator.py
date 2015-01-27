import os
import csv
import collections

website_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet/website/'


description_path = os.path.join(website_dir, 'download-descriptions.txt')
with open(description_path) as read_file:
    reader = csv.DictReader(read_file, delimiter='\t')
    rows = list(reader)

category_to_rows = collections.OrderedDict()

for row in rows:
    category_to_rows.setdefault(row['category'], list()).append(row)

row_html_template = (
'''
<div class="col-md-4">
  <h4>{name}</h4>
  <p>{description}</p>
  <p><a href="./files/{filename}" download>{filename}</a></p>
</div>
''')


html_list = list()
for category, category_rows in category_to_rows.items():
    html_str = '<h2>{}</h2>\n'.format(category)
    for i, row in enumerate(category_rows):
        if i % 3 == 0:
            html_str += '<div class="row">\n'
        html_str += row_html_template.format(**row)
        if i % 3 == 2:
            html_str += '</div>\n'
    if len(category_rows) % 3 != 0:
        html_str += '</div>\n'
    html_list.append(html_str)

html_str = '\n<hr>\n'.join(html_list)

description_path = os.path.join(website_dir, 'download-descriptions.html')
with open(description_path, 'w') as write_file:
    write_file.write(html_str)

# Description List Code
'''
#css
.dl-horizontal dt{
    width: 80px;
}
.dl-horizontal dd{
    margin-left: 90px;
}



#html
<dl class="dl-horizontal">
  <dt>Format</dt><dd>{format} </dd>
  <dt>Parsing</dt><dd>{parsing} </dd>
  <dt>Download</dt><dd><a href="./" download>{filename}</a></dd>
</dl>

'''