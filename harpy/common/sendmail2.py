#! /usr/bin/env python3

import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import socket

def send_pipeline_notification(status="completed", error_message=None):
    """
    Send an email notification when the pipeline finishes.
    Only sends if EMAIL environment variable is set.
    
    Args:
        status: Status of the pipeline (e.g., "completed", "failed")
        error_message: Optional error message if pipeline failed
    """
    email = os.environ.get('EMAIL')
    
    # Only proceed if EMAIL is set
    if not email:
        return
    
    # Extract domain from email for SMTP server lookup
    domain = email.split('@')[1] if '@' in email else None
    if not domain:
        print(f"Warning: Invalid email format in EMAIL variable: {email}")
        return
    
    # Construct email message
    msg = MIMEMultipart()
    msg['From'] = email
    msg['To'] = email
    msg['Subject'] = f"Pipeline {status.title()}"
    
    # Email body
    body = f"""
    Your pipeline has {status}.
    
    Hostname: {socket.gethostname()}
    Status: {status}
    """
    
    if error_message:
        body += f"\nError Details:\n{error_message}"
    
    msg.attach(MIMEText(body, 'plain'))
    
    # Try to send email using common SMTP ports
    smtp_servers = [
        (f'smtp.{domain}', 587),  # Standard submission port with STARTTLS
        (f'smtp.{domain}', 25),   # Standard SMTP port
        (f'mail.{domain}', 587),  # Alternative hostname
        (f'mail.{domain}', 25),   # Alternative hostname
    ]
    
    for server, port in smtp_servers:
        try:
            # Attempt connection
            with smtplib.SMTP(server, port, timeout=5) as smtp:
                smtp.ehlo()
                # Try STARTTLS if available
                if port == 587:
                    try:
                        smtp.starttls()
                        smtp.ehlo()
                    except:
                        pass
                
                # Send email (no authentication)
                smtp.send_message(msg)
                print(f"Pipeline notification sent to {email}")
                return
                
        except (smtplib.SMTPException, socket.error, OSError) as e:
            # Try next server/port combination
            continue
    
    # If all attempts failed, silently continue (don't break the pipeline)
    print(f"Note: Could not send email notification to {email}")


# Example usage:
if __name__ == "__main__":
    # At the end of your pipeline:
    try:
        # Your pipeline code here
        # ...
        send_pipeline_notification(status="completed")
    except Exception as e:
        send_pipeline_notification(status="failed", error_message=str(e))
        raise  # Re-raise the exception