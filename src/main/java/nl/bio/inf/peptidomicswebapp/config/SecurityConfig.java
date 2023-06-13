package nl.bio.inf.peptidomicswebapp.config;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.security.config.annotation.method.configuration.EnableMethodSecurity;
import org.springframework.security.config.annotation.web.builders.HttpSecurity;
import org.springframework.security.config.annotation.web.configuration.EnableWebSecurity;
import org.springframework.security.web.SecurityFilterChain;

/**
 *  This class builds the security filter chain for the application.
 * @author Jan Alfonso Busker
 */
@EnableWebSecurity
@EnableMethodSecurity
@Configuration
public class SecurityConfig {

    /**
     * Bean that configure the security chain
     * @param http HttpSecurity
     * @return HttpSecurity config
     * @throws Exception exception
     */
    @Bean
    protected SecurityFilterChain filterChain(HttpSecurity http) throws Exception {
        http
                .httpBasic().disable()
                .formLogin().disable()
                .csrf()
                .and()
                .cors();
        return http.build();
    }
}
