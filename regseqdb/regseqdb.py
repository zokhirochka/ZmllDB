import utils


############################################################################
# RegseqDB Class
class RegSeqDB:

	def __init__(self):
		'''
		RegseqDB Class Initializations
		'''

		# MariaDB Objects
		self.connection = None 
		self.cursor     = None 


	######################################
	# Database Methods
	def connect(self, host, port, database, username, password):
		"""Connect to Database and stores connection as class attribute.

		Args:
			host: Name of host database server
			port: Database server post 
			database: Database to connect
			uername: Username
			password: Password

		Returns:
			None
		"""
		keys = {
			"host": host,
			"port": port,
			"database": database,
			"username": username,
			"password": password
		}
		self.connection, self.cursor = utils.connect_db(keys)
		return

	def __query(self, query, inputs):
		"""Private generalized database query method

		Args:
			query: Structured SQL query template from calling function
			inputs: Curaated or cleaned inputs from calling functions

		Returns:
			Tuple of (results, colnames, rowcount)
		"""
		return utils.exec_query(
			cursor=self.cursor,
			query=query,
			inputs=inputs
		)


	######################################
	# Query Wrappers
	def get_promoter_expr(self, promoter, condition):
		"""Get Promoter Expression 

		Args:
			promoter: Name of promoter sequence to query
			condition: Condition in which to grab expression values

		Returns:
			Dictionary of {results, colnames, rowcount}
		"""

		# Check Promoter Value
		if not utils.db_contains(self.cursor, "Promoters", "pro_name", promoter):
			raise ValueError(f"ERROR: Promoter {promoter} not found in Promoters table.")
		# Check Condition Value
		if not utils.db_contains(self.cursor, "Experiments", "cond", condition):
			raise ValueError(f"ERROR: Condition {condition} not found in Experiments table.")

		# Define inputs and query
		inputs = [promoter, condition]
		query = """
				SELECT sID, num_DNA, num_RNA FROM Promoters
					JOIN PromoterSequences AS ps USING(pID)
					JOIN BarcodeCounts USING(sID)
					JOIN Experiments USING (eID)
				WHERE pID = (SELECT pID FROM Promoters WHERE pro_name = %s) AND 
					  eID = (SELECT eID FROM Experiments WHERE cond = %s LIMIT 1)
				GROUP BY ps.seq;
				"""

		# Submit Query
		results = self.__query(query, inputs)
		return results


	def get_promoter_expr_and_binding(self, promoter, condition, tf, include_rnap = True):
		"""Get Promoter Expression and TF/RNAP Binding

		Args:
			promoter: Name of promoter sequence to query
			condition: Condition in which to grab expression values
			tf: Name of TF to get binding
			include_rnap: Boolean to include RNAP binding

		Returns:
			Dictionary of {results, colnames, rowcount}
		"""

		# Check Promoter Value
		if not utils.db_contains(self.cursor, "Promoters", "pro_name", promoter):
			raise ValueError(f"ERROR: Promoter {promoter} not found in Promoters table.")
		# Check Condition Value
		if not utils.db_contains(self.cursor, "Experiments", "cond", condition):
			raise ValueError(f"ERROR: Condition {condition} not found in Experiments table.")
		# Check TF Value
		if not utils.db_contains(self.cursor, "TranscriptionFactors", "tf_name", tf):
			raise ValueError(f"ERROR: TF {tf} not found in TranscriptionFactors table.")

		# Define inputs and query
		inputs = [promoter, condition, tf]

		# Include RNAP Binding
		if include_rnap:
			query = """
					SELECT sID, SUM(num_DNA) as num_DNA, SUM(num_RNA) as num_RNA, energy, affinity FROM Promoters
						JOIN PromoterSequences AS ps USING(pID)
						JOIN BarcodeCounts USING(sID)
						JOIN Experiments USING (eID)
						JOIN BindingSitesRNAP using(sID)
						LEFT JOIN BindingSitesTF USING(sID)
						JOIN TranscriptionFactors USING(tID)
					WHERE pID = (SELECT pID FROM Promoters WHERE pro_name = %s)
						AND eID = (SELECT eID FROM Experiments WHERE cond = %s LIMIT 1)
						AND tID = (SELECT tID FROM TranscriptionFactors WHERE tf_name = %s)
					GROUP BY ps.seq;
					"""
		else:
			query = """
					SELECT sID, SUM(num_DNA) as num_DNA, SUM(num_RNA) as num_RNA, affinity FROM Promoters
						JOIN PromoterSequences AS ps USING(pID)
						JOIN BarcodeCounts USING(sID)
						JOIN Experiments USING (eID)
						JOIN BindingSitesTF USING(sID)
						JOIN TranscriptionFactors USING(tID)
					WHERE pID = (SELECT pID FROM Promoters WHERE pro_name = %s)
						AND eID = (SELECT eID FROM Experiments WHERE cond = %s LIMIT 1)
						AND tID = (SELECT tID FROM TranscriptionFactors WHERE tf_name = %s)
					GROUP BY ps.seq;
					"""

		# Submit Query
		results = self.__query(query, inputs)
		return results

