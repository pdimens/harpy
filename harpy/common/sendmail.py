#! /usr/bin/env python

import os
import sys
import smtplib
import subprocess
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

#def send_email(subject: str, message: str) -> None:
subject = "THIS"
message = "THAT"

try:
    you = os.environ["EMAIL"]
except KeyError:
    print("twas not found")
    sys.exit(1)
    #return

# Create the container (outer) email message.
me = 'harpy@cornell.edu'
subject = subject.strip('\"')
body = message.strip('\"')
attachment_type = 'plain'

if you.split('@')[0].isdigit():
    msg = MIMEText(body)
    msg['From'] = me
    msg['To'] = you
    msg['Subject'] = subject
else:
    msg = MIMEMultipart()
    msg.preamble = subject
    msg['From'] = me
    msg['To'] = you
    msg['Subject'] = subject
    msg.attach(MIMEText(body, attachment_type))

# Send the email via our own SMTP server.
try:
    s = smtplib.SMTP('localhost')
    s.sendmail(me, you, msg.as_string())
    s.quit()

except Exception as e:
    print("I FAILED:", e)
    sys.exit(1)
    #return
