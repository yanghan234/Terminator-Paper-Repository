#!/usr/bin/env python
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import sys

try:
    filename = sys.argv[1]
except IndexError:
    print("No Input File!!")
    exit()

fromaddr = "kaneyoung234@gmail.com"
toaddr = "kaneyoung234@gmail.com"

msg = MIMEMultipart()

msg['From'] = fromaddr
msg['To'] = toaddr
msg['Subject'] = "Fetch file "+filename+" from midway."

body = "The file is attached.\n"

msg.attach(MIMEText(body, 'plain'))

attachment = open(filename, "rb")

part = MIMEBase('application', 'octet-stream')
part.set_payload((attachment).read())
encoders.encode_base64(part)
part.add_header('Content-Disposition', "attachment; filename= %s" % filename)

msg.attach(part)

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login(fromaddr, "yanghan234")
text = msg.as_string()
server.sendmail(fromaddr, toaddr, text)
server.quit()

