
def log_handler(msg):
    if msg.get('level') == 'progress':
        cur = msg['done']
        total = msg['total']
        print(f"{cur} / {total} completed")
    elif ms.get('level') == 'job_info':
        curr = msg['msg']
        print(curr)
    