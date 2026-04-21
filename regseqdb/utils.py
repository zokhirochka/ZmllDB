import mariadb

############################################################################
# Database Methods

def connect_db(credentials_dictionary: dict) -> tuple:
        """Connect to MariaDB Database

        Args:
            credentials_dictionary: Dictionary containing required keys
                - host, port, databse, username, password

        Returns:
            mariadb.connect object

        Raises:
            ValueError: Missing required keys
        """

        # Validate Keys
        required_keys = ["host", "port", "database", "username", "password"]
        if not all(key in credentials_dictionary for key in required_keys):
                missing_keys = set(required_keys ) - set(credentials_dictionary.keys())
                raise ValueError(f"ERROR: Database Connection Error. Missing parameters {missing_keys}.")

        # Establish DB connection
        connection = mariadb.connect(**credentials_dictionary)
        cursor = connection.cursor()
        return (connection, cursor)


def exec_query(cursor, query, inputs):
        """Execute Query tonMariaDB Database

        Args:
            cursor: MariaDB cursor object
            query: Query Template
            input: Query Inputs

        Returns:
            Dictionary containing results, colnames, and rowcounts

        Raises:
            ValueError: Invalid Query
        """

        # Execute query and catch errors
        try:
            cursor.execute(query, inputs)
        except mariadb.Error as e:
            raise ValueError(f"ERROR: Invalid SQL Query.\nQuery: {query}\nError: {e}")

        # Fetch and Curate Results
        results_dict = {
                        "results": cursor.fetchall(),
                        "colnames": [metadata[0] for metadata in cursor.description],
                        "rowcount": len(results)
                        }
        return results_dict


def db_contains(cursor, table, column, value):
        """Checks if value is in MariaDB database

        Args:
            cursor: MariaDB cursor object
            table: Table to serch for value
            column: Column name for value search
            value: Value to match in column

        Returns:
            True or False based on value presence
        """
        inputs = [value]
        query = f"SELECT * FROM {table} WHERE {column} = %s;"
        results = exec_query(cursor=cursor, query=query, inputs=inputs)
        if not results:
                return False
        return True


