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


results = db.get_promoter_expr_and_binding("tisBp", "stationary phase (3d)", "LexA", include_rnap=True)
print("\t".join(results["colnames"]))
for r in results["results"]:
        tmp = []
        for t in r:
                tmp.append(str(t))
        print("\t".join(tmp))