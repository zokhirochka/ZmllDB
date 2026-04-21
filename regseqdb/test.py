from regseqdb import RegSeqDB

username = str(input("Provide username: "))
password = str(input("Provide password: "))

db = RegSeqDB()
db.connect(
	host="bioed-new.bu.edu",
	port=4253,
	database="Team12",
	username=username,
	password=password
)
